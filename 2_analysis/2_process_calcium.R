# =============================================================
# 2_process_calcium.R
# Process per-cell calcium transient data from CyteSeer traces
#
# Loads its OWN data from 1_data/Calcium/ (2 dirs: 1.6, 1.6.24)
# Uses shared helpers from 0_setup/helpers.R:
#   load_plate_maps, load_and_merge_cyteseer
# Extracts calcium params + Cell.Area + Nucleus.Area + NC_ratio
#
# Expects in global env:
#   - main_dir, analysis_date, group_levels
#   - Packages already loaded (tidyverse, pracma, etc.)
#   - Utility functions from contraction_trace_utils.R
# =============================================================

cat("=== 2_process_calcium.R ===\n")

# ---- Output directories ----
qc_dir        <- file.path(main_dir, "QC")
processed_dir <- file.path(main_dir, "processed")
dir.create(qc_dir,        showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Pixel size for area conversion ----
# 0.325 µm/pixel -> area in µm² = (pixels²) * 0.325²
px_size_um <- 0.325

# ---- 1. Load per-cell data from CyteSeer (Calcium session) ----
cat("Loading CyteSeer single-cell data (Calcium session)...\n")

cyteseer_dirs <- tibble::tribble(
  ~plate_id,  ~root_dir,
  "1.6",      "../1_data/Calcium/CyteSeerResults_1.6_B02_B04",
  "1.6",      "../1_data/Calcium/CyteSeerResults_1.6_B05_C04",
  "1.6",      "../1_data/Calcium/CyteSeerResults_1.6_C05_E02",
  "1.6",      "../1_data/Calcium/CyteSeerResults_1.6_E03_E07",
  "1.6.24",   "../1_data/Calcium/CyteSeerResults_1.6.24_B02_B11",
  "1.6.24",   "../1_data/Calcium/CyteSeerResults_1.6.24_C02_C11",
  "1.6.24",   "../1_data/Calcium/CyteSeerResults_1.6.24_D02_D06",
  "1.6.24",   "../1_data/Calcium/CyteSeerResults_1.6.24_D07_D11",
  "1.6.24",   "../1_data/Calcium/CyteSeerResults_1.6.24_E02_E11"
)

plate_map_all <- load_plate_maps()
singlecell <- load_and_merge_cyteseer(cyteseer_dirs, plate_map_all, skip = 224)
cat(sprintf("  Starting with %d cells\n", nrow(singlecell)))

# ---- 2. Parse and analyse calcium traces ----
# Use the CyteSeer "Trimmed" trace (already a normalised ΔF/F-style
# signal) with NO further preprocessing — no smoothing, no detrending,
# no baseline subtraction. See analyze_cell_calcium().
cat("  Processing calcium traces (Trimmed, no preproc)...\n")
ca_results <- lapply(singlecell$Trimmed.Ci.on.Cm.API.F.t.,
                     analyze_cell_calcium)
cat(sprintf("  Processed %d calcium traces\n", length(ca_results)))

# Count how many returned valid results
n_valid <- sum(sapply(ca_results, function(x) !is.null(x) && !is.null(x$params)))
n_null  <- sum(sapply(ca_results, is.null))
cat(sprintf("  Valid params: %d | NULL trace: %d | No params: %d\n",
            n_valid, n_null, length(ca_results) - n_valid - n_null))

# ---- 3. Collect per-cell calcium parameters ----
cat("  Collecting per-cell calcium parameters...\n")

valid_idx <- which(sapply(ca_results, function(x) !is.null(x) && !is.null(x$params)))
if (length(valid_idx) == 0) stop("No valid calcium parameters found - check traces")

ca_param_names <- names(ca_results[[valid_idx[1]]]$params)

# CyteSeer's own per-cell "is Responder" flag (does this cell have a detectable
# transient at all). Carried through so run_pipeline.R can gate wells on the
# fraction of responding cells — a signal-quality measure computed by the
# instrument software, independent of the trace analysis in this script.
resp_col <- "Ci.on.Cm.API.is.Responder"
if (!resp_col %in% names(singlecell)) {
  stop("Column '", resp_col, "' not found in the CyteSeer export. It is needed ",
       "for the well-level responder QC gate in run_pipeline.R.")
}
responder_flag <- tolower(trimws(as.character(singlecell[[resp_col]]))) %in% "true"

ca_cell_params <- do.call(rbind, lapply(seq_along(ca_results), function(i) {
  res <- ca_results[[i]]
  if (is.null(res) || is.null(res$params)) {
    params <- as.list(rep(NA_real_, length(ca_param_names)))
    names(params) <- ca_param_names
  } else {
    params <- res$params
  }
  data.frame(
    cell_uid     = singlecell$cell_uid[i],
    well_uid     = singlecell$well_uid[i],
    Plate        = singlecell$Plate[i],
    Well         = singlecell$Well[i],
    Condition    = singlecell$Condition[i],
    Cell_Line    = singlecell$Cell_Line[i],
    Cell.Area    = singlecell$Area.Cm[i] * px_size_um^2,
    Nucleus.Area = singlecell$Area.Nm[i] * px_size_um^2,
    NC_ratio     = singlecell$Area.Nm[i] / (singlecell$Area.Cm[i] + 1e-6),
    is_responder = responder_flag[i],
    params,
    row.names = NULL
  )
}))

cat(sprintf("  Calcium parameter table: %d rows x %d cols\n",
            nrow(ca_cell_params), ncol(ca_cell_params)))

# ---- 4. QC classification — AC1 pacing score + hard n_beats cap ----
# Single data-driven rule:
#   good_paced = (AC1 >= Otsu(AC1)) AND (n_beats <= max_beats_ca)
# AC1 is the normalised autocorrelation at lag = one pacing period (1 s).
# max_beats_ca = 10 is the number of beats delivered by the 1 Hz pacemaker
# over the 10 s calcium recording. Any cell with more beats is by definition
# not following pacing (spontaneous / ectopic). Fewer beats is allowed (may
# reflect a missed peak due to bleaching). Everything that fails either test
# is classified "irregular". See 0_setup/contraction_trace_utils.R for
# score_pacing_ac1() and otsu_threshold().

cat("  Applying QC (AC1 + n_beats cap)...\n")

max_beats_ca <- 10   # 1 Hz pacing, 10 s recording

# AC1 per cell from the trace list (aligned with ca_cell_params row order)
ca_cell_params$ac1 <- vapply(
  ca_results,
  function(r) if (is.null(r)) NA_real_ else score_pacing_ac1(r$trace),
  numeric(1)
)

# Data-driven threshold on the AC1 distribution
ac1_otsu_ca <- otsu_threshold(ca_cell_params$ac1)
cat(sprintf("  AC1 Otsu threshold: %.4f\n", ac1_otsu_ca))

ca_cell_params_qc <- ca_cell_params %>%
  mutate(
    qc_pass = !is.na(n_beats),
    signal_class = case_when(
      !qc_pass                                          ~ "irregular",
      is.finite(ac1) & ac1 >= ac1_otsu_ca &
        n_beats <= max_beats_ca                         ~ "good_paced",
      TRUE                                              ~ "irregular"
    )
  )

## Print summary
n_total    <- nrow(ca_cell_params_qc)
n_good     <- sum(ca_cell_params_qc$signal_class == "good_paced")
n_irreg    <- n_total - n_good
n_cap_drop <- sum(is.finite(ca_cell_params_qc$ac1) &
                    ca_cell_params_qc$ac1 >= ac1_otsu_ca &
                    ca_cell_params_qc$n_beats > max_beats_ca,
                  na.rm = TRUE)
cat(sprintf("  good_paced: %d / %d (%.1f%%)\n",
            n_good, n_total, 100 * n_good / n_total))
cat(sprintf("  irregular : %d / %d (%.1f%%)\n",
            n_irreg, n_total, 100 * n_irreg / n_total))
cat(sprintf("    of which dropped ONLY by n_beats > %d cap: %d\n",
            max_beats_ca, n_cap_drop))

# ---- 5. QC visualizations ----
cat("  Generating calcium QC plots...\n")

pdf_path <- file.path(qc_dir, "calcium_qc_plots.pdf")
cairo_pdf(pdf_path, width = 16, height = 10, onefile = TRUE)

# --- 5a. Per-well calcium trace overlay (MOST IMPORTANT plot) ---
cat("    Per-well trace overlay...\n")

# Build a data frame of all traces with well info
trace_list <- lapply(seq_along(ca_results), function(i) {
  res <- ca_results[[i]]
  if (is.null(res) || is.null(res$trace)) return(NULL)
  res$trace %>%
    mutate(
      cell_uid = singlecell$cell_uid[i],
      well_uid = singlecell$well_uid[i],
      Condition = singlecell$Condition[i],
      Cell_Line = singlecell$Cell_Line[i]
    )
})
trace_df <- bind_rows(trace_list)

if (nrow(trace_df) > 0) {
  # Get unique wells and plot in pages of ~20 wells each
  wells <- sort(unique(trace_df$well_uid))
  n_per_page <- 20
  n_pages <- ceiling(length(wells) / n_per_page)

  for (pg in seq_len(n_pages)) {
    idx_start <- (pg - 1) * n_per_page + 1
    idx_end   <- min(pg * n_per_page, length(wells))
    wells_page <- wells[idx_start:idx_end]

    p_traces <- trace_df %>%
      filter(well_uid %in% wells_page) %>%
      ggplot(aes(x = time_ms, y = signal_bc, group = cell_uid)) +
      geom_line(alpha = 0.15, linewidth = 0.2) +
      facet_wrap(~ well_uid, scales = "free_y") +
      labs(
        title = sprintf("Calcium traces per well (page %d/%d)", pg, n_pages),
        x = "Time (ms)",
        y = "Calcium signal (CyteSeer Trimmed)"
      ) +
      theme_project()

    print(p_traces)
  }
}

# --- 5b. AC1 distribution with Otsu threshold ---
cat("    Plotting AC1 distribution...\n")

class_colors <- c(good_paced = "#0072B2", irregular = "#D55E00")

p_ac1 <- ca_cell_params_qc %>%
  filter(is.finite(ac1)) %>%
  ggplot(aes(x = ac1, fill = signal_class)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
  geom_vline(xintercept = ac1_otsu_ca, linetype = "dashed",
             linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title = sprintf(
         "Calcium AC1 distribution (Otsu threshold = %.3f)",
         ac1_otsu_ca),
       subtitle = sprintf("good_paced also requires n_beats <= %d",
                          max_beats_ca),
       x = "AC1 (autocorrelation @ 1 s lag)", y = "Cells",
       fill = NULL) +
  theme_project() +
  theme(legend.position = "bottom")
print(p_ac1)

# --- 5c. n_beats distribution with the max_beats cap ---
p_nbeats <- ca_cell_params_qc %>%
  filter(!is.na(n_beats)) %>%
  ggplot(aes(x = n_beats, fill = signal_class)) +
  geom_histogram(binwidth = 1, colour = "white", alpha = 0.9) +
  geom_vline(xintercept = max_beats_ca + 0.5,
             linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title = sprintf(
         "Calcium n_beats (cap = %d, cells with > cap are irregular)",
         max_beats_ca),
       x = "n_beats", y = "Cells", fill = NULL) +
  theme_project() +
  theme(legend.position = "bottom")
print(p_nbeats)

# --- 5d. AC1 vs n_beats scatter (QC decision plane) ---
p_scatter <- ca_cell_params_qc %>%
  filter(is.finite(ac1), !is.na(n_beats)) %>%
  ggplot(aes(x = ac1, y = n_beats, colour = signal_class)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = ac1_otsu_ca,
             linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = max_beats_ca + 0.5,
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller) +
  scale_colour_manual(values = class_colors) +
  labs(title = "Calcium: AC1 vs n_beats (QC decision plane)",
       x = "AC1", y = "n_beats", colour = NULL) +
  theme_project() +
  theme(legend.position = "bottom")
print(p_scatter)

# --- 5e. Example traces (good_paced + irregular) ---
cat("    Example traces...\n")
set.seed(42)

sample_examples <- function(cls, n = 6) {
  idx <- which(ca_cell_params_qc$signal_class == cls)
  if (length(idx) == 0) return(tibble())
  if (length(idx) > n) idx <- sample(idx, n)
  bind_rows(lapply(idx, function(i) {
    res <- ca_results[[i]]
    if (is.null(res) || is.null(res$trace)) return(NULL)
    ci <- ca_cell_params_qc[i, ]
    label <- sprintf(
      "%s\n%s | %s\nn=%s  AC1=%.3f  medIBI=%.0f",
      ci$cell_uid, ci$Condition, cls,
      ifelse(is.na(ci$n_beats), "NA", as.character(ci$n_beats)),
      ifelse(is.na(ci$ac1), NA, ci$ac1),
      ifelse(is.na(ci$median_ibi_ms), NA, ci$median_ibi_ms)
    )
    res$trace %>%
      mutate(facet_label = label, example_group = cls)
  }))
}

example_trace_df <- bind_rows(
  sample_examples("good_paced"),
  sample_examples("irregular")
)

if (nrow(example_trace_df) > 0) {
  for (grp in unique(example_trace_df$example_group)) {
    grp_df <- example_trace_df %>% filter(example_group == grp)
    if (nrow(grp_df) == 0) next
    p_ex <- ggplot(grp_df, aes(x = time_ms, y = signal_bc)) +
      geom_line(linewidth = 0.3) +
      facet_wrap(~ facet_label, scales = "free_y") +
      labs(title = sprintf("Calcium: example traces — %s", grp),
           x = "Time (ms)", y = "Signal") +
      theme_project() +
      theme(strip.text = element_text(size = 7))
    print(p_ex)
  }
}

dev.off()
cat(sprintf("  Saved QC plots: %s\n", pdf_path))

# ---- Well-level QC summary table ----
# Reports the real QC outcome: fraction of cells classified as "good_paced"
# (AC1 >= Otsu + n_beats <= cap). `qc_pass` alone only marks whether the cell
# had detectable beats and is near-100% by construction, so it is not used
# here for the pass rate.
ca_well_qc <- ca_cell_params_qc %>%
  group_by(Plate, Well, well_uid, Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = sum(signal_class == "irregular"),
    pct_pass     = round(100 * n_good_paced / n_total, 1),
    .groups  = "drop"
  ) %>%
  arrange(pct_pass)

well_qc_path <- file.path(qc_dir, "calcium_well_qc_summary.csv")
write.csv(ca_well_qc, well_qc_path, row.names = FALSE)
cat(sprintf("  Saved well QC summary: %s\n", well_qc_path))

# ---- 5b. QC Calibration gallery (stratified ~50 cells for labeling) ----

cat("  Generating QC calibration gallery...\n")
set.seed(42)

# Stratify around the AC1 Otsu cutoff so the gallery focuses on the
# decision boundary (cells with AC1 just above / just below the threshold).
gp_cells  <- ca_cell_params_qc %>% filter(signal_class == "good_paced")
irr_cells <- ca_cell_params_qc %>% filter(signal_class == "irregular")

strat1 <- gp_cells %>%
  filter(is.finite(ac1)) %>%
  slice_min(ac1, n = 15, with_ties = FALSE) %>%
  mutate(stratum = "gp_just_above_otsu")
strat2 <- irr_cells %>%
  filter(is.finite(ac1)) %>%
  slice_max(ac1, n = 10, with_ties = FALSE) %>%
  mutate(stratum = "irr_just_below_otsu")
strat3 <- irr_cells %>%
  filter(is.finite(ac1),
         !is.na(n_beats) & n_beats > max_beats_ca) %>%
  slice_sample(n = min(10, nrow(.))) %>%
  mutate(stratum = "irr_nbeats_cap")
strat4 <- gp_cells %>%
  filter(!cell_uid %in% strat1$cell_uid) %>%
  slice_sample(n = min(10, nrow(.))) %>%
  mutate(stratum = "gp_random")
strat5 <- irr_cells %>%
  filter(!cell_uid %in% c(strat2$cell_uid, strat3$cell_uid)) %>%
  slice_sample(n = min(5, nrow(.))) %>%
  mutate(stratum = "irr_random")

calib_cells <- bind_rows(strat1, strat2, strat3, strat4, strat5) %>%
  distinct(cell_uid, .keep_all = TRUE) %>%
  mutate(id = row_number())

n_calib <- nrow(calib_cells)
cat("  Calibration sample:", n_calib, "cells\n")

# Save companion CSV (calcium-specific: include mean_amp, CTD50, CTD90)
calib_csv <- calib_cells %>%
  select(id, cell_uid, well_uid, Condition, Cell_Line,
         signal_class, stratum, ac1, n_beats, median_ibi_ms,
         snr_amp, rise_time_ms, decay_time_ms,
         mean_amp, CTD50_ms, CTD90_ms) %>%
  mutate(user_label = "")
calib_csv_path <- file.path(qc_dir,
  paste0("calibration_calcium_", analysis_date, ".csv"))
write.csv(calib_csv, calib_csv_path, row.names = FALSE)
cat("  Calibration CSV saved to:", calib_csv_path, "\n")

# Build trace data with peak markers
calib_trace_list <- lapply(seq_len(n_calib), function(j) {
  uid     <- calib_cells$cell_uid[j]
  row_idx <- which(singlecell$cell_uid == uid)[1]
  if (is.na(row_idx)) return(NULL)
  res <- ca_results[[row_idx]]
  if (is.null(res) || is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
  ci <- calib_cells[j, ]
  label <- sprintf(
    "#%d | %s | %s | %s\nn=%s  AC1=%.3f  IBI=%.0f  snr=%.1f",
    ci$id, ci$cell_uid, ci$signal_class, ci$stratum,
    ifelse(is.na(ci$n_beats), "NA", as.character(ci$n_beats)),
    ifelse(is.na(ci$ac1), NA, ci$ac1),
    ifelse(is.na(ci$median_ibi_ms), NA, ci$median_ibi_ms),
    ifelse(is.na(ci$snr_amp), NA, ci$snr_amp)
  )
  trace_df <- tibble(
    id           = ci$id,
    facet_label  = label,
    time_ms      = res$trace$time_ms,
    signal_bc    = res$trace$signal_bc,
    is_peak      = FALSE
  )
  if (!is.null(res$peaks) && nrow(res$peaks) > 0) {
    peak_df <- tibble(
      id           = ci$id,
      facet_label  = label,
      time_ms      = res$peaks$peak_time,
      signal_bc    = res$peaks$peak_amp,
      is_peak      = TRUE
    )
    trace_df <- bind_rows(trace_df, peak_df)
  }
  trace_df
})
calib_trace_df <- bind_rows(calib_trace_list)

if (nrow(calib_trace_df) > 0) {
  calib_trace_df$facet_label <- factor(
    calib_trace_df$facet_label,
    levels = unique(calib_trace_df$facet_label[
      order(calib_trace_df$id)]))

  traces_only <- calib_trace_df %>% filter(!is_peak)
  peaks_only  <- calib_trace_df %>% filter(is_peak)

  calib_pdf_path <- file.path(qc_dir,
    paste0("calibration_calcium_", analysis_date, ".pdf"))
  n_per_page <- 6
  page_ids   <- split(sort(unique(calib_cells$id)),
    ceiling(seq_along(unique(calib_cells$id)) / n_per_page))

  cairo_pdf(calib_pdf_path, width = 14, height = 10, onefile = TRUE)
  for (pg in page_ids) {
    pg_labels <- levels(calib_trace_df$facet_label)[pg]
    pg_traces <- traces_only %>%
      filter(facet_label %in% pg_labels)
    pg_peaks  <- peaks_only %>%
      filter(facet_label %in% pg_labels)

    p <- ggplot(pg_traces, aes(x = time_ms, y = signal_bc)) +
      geom_line(linewidth = 0.4, color = "grey30") +
      geom_point(data = pg_peaks, color = "red", size = 2) +
      facet_wrap(~ facet_label, scales = "free",
                 ncol = 3, nrow = 2) +
      labs(x = "Time (ms)",
           y = "Calcium signal (CyteSeer Trimmed)") +
      theme_project() +
      theme(strip.text = element_text(size = 7))
    print(p)
  }
  dev.off()
  cat("  Calibration PDF saved to:", calib_pdf_path, "\n")
}

# ---- 5c. Publication trace gallery ----

cat("  Generating trace gallery...\n")

if (!exists("expl_dir"))
  expl_dir <- file.path(main_dir, "Figures", "Exploratory")
gallery_dir <- file.path(expl_dir, "trace_galleries")
dir.create(gallery_dir, showWarnings = FALSE, recursive = TRUE)

## Pick 5 representative wells (top 3 per cell line, then trim to 5) ----------
gallery_well_pick <- ca_cell_params_qc %>%
  filter(signal_class == "good_paced") %>%
  count(well_uid, Cell_Line, name = "n_good_paced") %>%
  group_by(Cell_Line) %>%
  slice_max(n_good_paced, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  slice_max(n_good_paced, n = 5, with_ties = FALSE)

set.seed(123)

gallery_colors <- c("good_paced" = "#0072B2",
                    "irregular"  = "#D55E00")

## One PDF per well, 16 traces each (4 cols x 4 rows) ----------------------
for (w in seq_len(nrow(gallery_well_pick))) {
  wuid <- gallery_well_pick$well_uid[w]
  pool <- ca_cell_params_qc %>% filter(well_uid == wuid)

  gp  <- pool %>% filter(signal_class == "good_paced")
  irr <- pool %>% filter(signal_class == "irregular")

  picked <- bind_rows(
    if (nrow(gp)  > 0) slice_sample(gp,  n = min(12, nrow(gp)))  else NULL,
    if (nrow(irr) > 0) slice_sample(irr, n = min(4,  nrow(irr))) else NULL
  )
  if (nrow(picked) > 16) picked <- slice_sample(picked, n = 16)

  trace_list <- lapply(seq_len(nrow(picked)), function(j) {
    uid     <- picked$cell_uid[j]
    row_idx <- which(singlecell$cell_uid == uid)[1]
    if (is.na(row_idx)) return(NULL)
    res <- ca_results[[row_idx]]
    if (is.null(res) || is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
    ci    <- picked[j, ]
    label <- sprintf(
      "%s\n%s | n=%s  AC1=%.3f  medIBI=%.0f  snr=%.1f",
      uid, ci$signal_class,
      ifelse(is.na(ci$n_beats), "NA", as.character(ci$n_beats)),
      ifelse(is.na(ci$ac1), NA, ci$ac1),
      ifelse(is.na(ci$median_ibi_ms), NA, ci$median_ibi_ms),
      ifelse(is.na(ci$snr_amp), NA, ci$snr_amp)
    )
    tibble(facet_label  = label,
           signal_class = ci$signal_class,
           time_ms      = res$trace$time_ms,
           signal_bc    = res$trace$signal_bc)
  })
  trace_df <- bind_rows(trace_list)

  if (nrow(trace_df) == 0) next

  p <- ggplot(trace_df,
              aes(x = time_ms, y = signal_bc, color = signal_class)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ facet_label, scales = "free", ncol = 4) +
    scale_color_manual(values = gallery_colors) +
    labs(title    = paste("Calcium gallery:", wuid),
         x = "Time (ms)", y = "Calcium signal (CyteSeer Trimmed)",
         color = "Signal class") +
    theme_project() +
    theme(strip.text      = element_text(size = 6),
          legend.position = "bottom")

  safe_wuid <- gsub("[^A-Za-z0-9_.]", "_", wuid)
  pdf_path  <- file.path(
    gallery_dir,
    paste0("trace_gallery_calcium_", safe_wuid,
           "_", analysis_date, ".pdf")
  )
  cairo_pdf(pdf_path, width = 14, height = 12)
  print(p)
  dev.off()
  cat("  Gallery saved:", pdf_path, "\n")
}

# ---- 6. Save outputs ----
rds_path <- file.path(processed_dir, "calcium_all_cells.rds")
saveRDS(ca_cell_params_qc, rds_path)
cat(sprintf("  Saved calcium RDS: %s\n", rds_path))

## RDS with raw trace results (needed for cluster trace plots when skipping step 2)
results_rds_path <- file.path(processed_dir, "calcium_results_all.rds")
saveRDS(ca_results, results_rds_path)
cat(sprintf("  Saved calcium results RDS: %s\n", results_rds_path))

# Store in global env for downstream scripts
ca_cell_params_all  <<- ca_cell_params_qc
ca_results_all      <<- ca_results
ca_well_qc_summary  <<- ca_well_qc

cat(sprintf("=== 2_process_calcium.R complete | %d cells, %d PASS (%.1f%%) ===\n",
            nrow(ca_cell_params_qc),
            sum(ca_cell_params_qc$qc_pass),
            100 * mean(ca_cell_params_qc$qc_pass)))
