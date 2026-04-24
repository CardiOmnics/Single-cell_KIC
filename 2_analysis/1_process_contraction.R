###############################################################################
#  1_process_contraction.R
#  -------------------------------------------------------------------------
#  Process per-cell contraction data from CyteSeer Primary_DataTable.csv files.
#
#  Uses shared helpers from 0_setup/helpers.R:
#    - read_cyteseer_singlecell()  (read a single CyteSeer folder)
#    - load_plate_maps()           (load and combine plate maps)
#    - load_and_merge_cyteseer()   (read + merge + add identifiers)
#
#  Loads from 1_data/Contraction/ (3 dirs: 1.6, 1.6.24_B, 1.6.24_CDE).
#  Extracts contraction parameters + Bodipy intensity + Nucleus.Area.
#
#  Sourced from the master pipeline script which has already loaded:
#    - 0_setup/colors_1.R   (group_colors, group_levels)
#    - 0_setup/theme.R      (theme_project())
#    - 0_setup/helpers.R
#    - 0_setup/contraction_trace_utils.R
#        (parse_trace_string, preprocess_trace, detect_beats,
#         extract_paced_contraction_parameters, analyze_well_trace,
#         extract_beat_kinetics)
#
#  Master script also provides: main_dir, analysis_date
#  Packages already loaded: tidyverse, ggpubr, broom, lme4, emmeans,
#                            openxlsx, lmerTest, broom.mixed, pracma
###############################################################################

cat("=== 1_process_contraction.R ===\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load per-cell data from CyteSeer folders
# ─────────────────────────────────────────────────────────────────────────────

cat("Loading CyteSeer single-cell data...\n")

## Pixel size for area conversion -------------------------------------------
# 0.325 µm/pixel -> area in µm² = (pixels²) * 0.325²
px_size_um <- 0.325

## Directories containing CyteSeer output ----------------------------------
cyteseer_dirs <- tibble::tribble(
  ~plate_id,  ~root_dir,
  "1.6",      "../1_data/Contraction/CyteSeerResults_2025-11-17_15-38-19_1.6",
  "1.6.24",   "../1_data/Contraction/CyteSeerResults_2025-11-17_16-56-10_1.6.24_B",
  "1.6.24",   "../1_data/Contraction/CyteSeerResults_2025-12-04_12-50-13_1.6.24_CDE"
)

## Load plate maps and single-cell data using shared helpers ----------------
plate_map_all <- load_plate_maps()
singlecell <- load_and_merge_cyteseer(cyteseer_dirs, plate_map_all)

# ─────────────────────────────────────────────────────────────────────────────
# 1b. Extract Bodipy intensity from TMRM channel
# ─────────────────────────────────────────────────────────────────────────────
# Bodipy was stained in the contraction session. TMRM channel maps to Bodipy.
# It's a static stain so the time series is repeated measurements — take mean.
cat("  Extracting Bodipy intensity from TMRM channel...\n")
bodipy_intensity <- sapply(singlecell$TMRM.API.F.t., function(x) {
  df <- parse_trace_string(x)
  if (is.null(df)) return(NA_real_)
  mean(df$signal, na.rm = TRUE)
})
singlecell$Bodipy_intensity <- bodipy_intensity
cat(sprintf("    Bodipy values: %d of %d cells\n",
            sum(!is.na(bodipy_intensity)), length(bodipy_intensity)))

# ─────────────────────────────────────────────────────────────────────────────
# 2. Parse contraction traces
# ─────────────────────────────────────────────────────────────────────────────

cat("Processing contraction traces...\n")
# normalized = TRUE -> skip preprocess_trace() entirely: no smoothing,
# no linear detrend, no baseline subtraction. signal_bc <- signal as-is.
contr_results <- lapply(singlecell$Calcium.Stretch.F.t.,
                        analyze_well_trace, normalized = TRUE)
cat("  Traces processed:", length(contr_results), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Collect per-cell parameters
# ─────────────────────────────────────────────────────────────────────────────

cat("Collecting per-cell contraction parameters...\n")

valid_idx <- which(sapply(contr_results, function(x) !is.null(x$params)))
if (length(valid_idx) == 0) stop("No valid contraction params found")
param_names <- names(contr_results[[valid_idx[1]]]$params)

contr_cell_params <- do.call(rbind, lapply(seq_along(contr_results), function(i) {
  res <- contr_results[[i]]
  if (is.null(res$params)) {
    params <- as.list(rep(NA_real_, length(param_names)))
    names(params) <- param_names
  } else {
    params <- res$params
  }
  data.frame(
    cell_uid            = singlecell$cell_uid[i],
    well_uid            = singlecell$well_uid[i],
    Plate               = singlecell$Plate[i],
    Well                = singlecell$Well[i],
    Condition           = singlecell$Condition[i],
    Cell_Line           = singlecell$Cell_Line[i],
    Nucleus.Area        = singlecell$Nucleus.Area[i] * px_size_um^2,
    Bodipy_intensity    = singlecell$Bodipy_intensity[i],
    params,
    row.names = NULL
  )
}))

cat("  Parameter table rows:", nrow(contr_cell_params), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. QC classification — AC1 pacing score + hard n_beats cap
# ─────────────────────────────────────────────────────────────────────────────
# Single data-driven rule:
#   good_paced = (AC1 >= Otsu(AC1)) AND (n_beats <= max_beats_con)
# AC1 is the normalised autocorrelation at lag = one pacing period (1 s).
# max_beats_con = 5 is the number of beats delivered by the 1 Hz pacemaker
# over the ~5 s contraction recording — any cell showing more beats is by
# definition not following pacing (spontaneous / ectopic). Fewer beats is
# allowed (may reflect a missed detection). Everything that fails either
# test is classified "irregular" (there is no "poor_signal" class; noisy
# traces are biologically indistinguishable from arrhythmic cells here).
# See 0_setup/contraction_trace_utils.R for score_pacing_ac1() and
# otsu_threshold().

cat("Applying QC (AC1 + n_beats cap)...\n")

max_beats_con <- 5     # 1 Hz pacing, ~5 s recording

# AC1 per cell from the trace list (aligned with contr_cell_params row order)
contr_cell_params$ac1 <- vapply(
  contr_results,
  function(r) score_pacing_ac1(r$trace),
  numeric(1)
)

# Data-driven threshold on the AC1 distribution
ac1_otsu_con <- otsu_threshold(contr_cell_params$ac1)
cat(sprintf("  AC1 Otsu threshold: %.4f\n", ac1_otsu_con))

contr_cell_params_qc <- contr_cell_params %>%
  mutate(
    # qc_pass = has valid data to classify (keeps the cell in the table).
    # Downstream code uses filter(qc_pass) to drop rows with missing params
    # before further splitting by signal_class.
    qc_pass = !is.na(n_beats),
    signal_class = case_when(
      !qc_pass                                         ~ "irregular",
      is.finite(ac1) & ac1 >= ac1_otsu_con &
        n_beats <= max_beats_con                       ~ "good_paced",
      TRUE                                             ~ "irregular"
    )
  )

## Print summary
n_total   <- nrow(contr_cell_params_qc)
n_good    <- sum(contr_cell_params_qc$signal_class == "good_paced")
n_irreg   <- n_total - n_good
n_cap_drop <- sum(is.finite(contr_cell_params_qc$ac1) &
                    contr_cell_params_qc$ac1 >= ac1_otsu_con &
                    contr_cell_params_qc$n_beats > max_beats_con,
                  na.rm = TRUE)
cat(sprintf("  good_paced: %d / %d (%.1f%%)\n",
            n_good, n_total, 100 * n_good / n_total))
cat(sprintf("  irregular : %d / %d (%.1f%%)\n",
            n_irreg, n_total, 100 * n_irreg / n_total))
cat(sprintf("    of which dropped ONLY by n_beats > %d cap: %d\n",
            max_beats_con, n_cap_drop))

# ─────────────────────────────────────────────────────────────────────────────
# 5. QC visualizations - multi-page PDF
# ─────────────────────────────────────────────────────────────────────────────

cat("Generating QC plots...\n")

## Create output directories -----------------------------------------------
qc_dir        <- file.path(main_dir, "QC")
processed_dir <- file.path(main_dir, "processed")
dir.create(qc_dir,        showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir,  showWarnings = FALSE, recursive = TRUE)

## Open multi-page PDF ------------------------------------------------------
pdf_path <- file.path(qc_dir, "contraction_qc_plots.pdf")
cairo_pdf(pdf_path, width = 14, height = 10, onefile = TRUE)

# --- Plot 1: Per-well trace overlay (key diagnostic) ----------------------
cat("  Building per-well trace overlay...\n")

# Build long-form trace table from contr_results
trace_list <- lapply(seq_along(contr_results), function(i) {
  res <- contr_results[[i]]
  if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
  tibble(
    cell_uid  = singlecell$cell_uid[i],
    well_uid  = singlecell$well_uid[i],
    Condition = as.character(singlecell$Condition[i]),
    Cell_Line = as.character(singlecell$Cell_Line[i]),
    time_ms   = res$trace$time_ms,
    signal_bc = res$trace$signal_bc
  )
})
trace_long <- bind_rows(trace_list)

if (nrow(trace_long) > 0) {
  wells_all <- sort(unique(trace_long$well_uid))
  wells_per_page <- 12
  n_pages <- ceiling(length(wells_all) / wells_per_page)

  for (pg in seq_len(n_pages)) {
    idx_start <- (pg - 1) * wells_per_page + 1
    idx_end   <- min(pg * wells_per_page, length(wells_all))
    wells_sub <- wells_all[idx_start:idx_end]

    p_traces <- trace_long %>%
      filter(well_uid %in% wells_sub) %>%
      ggplot(aes(x = time_ms, y = signal_bc, group = cell_uid)) +
      geom_line(alpha = 0.15, linewidth = 0.2) +
      facet_wrap(~ well_uid, scales = "free_y", ncol = 4) +
      labs(title = paste0("Contraction trace overlay (page ", pg, "/", n_pages, ")"),
           x = "Time (ms)", y = "Contraction (raw stretch)") +
      theme_project() +
      theme(strip.text = element_text(size = 7))
    print(p_traces)
  }
} else {
  plot.new()
  text(0.5, 0.5, "No contraction traces available for overlay plot", cex = 1.2)
}

# --- Plot 2: AC1 distribution with Otsu threshold --------------------------
cat("  Plotting AC1 distribution...\n")

class_colors <- c(good_paced = "#2166AC", irregular = "#D95F02")

p_ac1 <- contr_cell_params_qc %>%
  filter(is.finite(ac1)) %>%
  ggplot(aes(x = ac1, fill = signal_class)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
  geom_vline(xintercept = ac1_otsu_con, linetype = "dashed",
             linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title = sprintf(
         "Contraction AC1 distribution (Otsu threshold = %.3f)",
         ac1_otsu_con),
       subtitle = sprintf("good_paced also requires n_beats <= %d",
                          max_beats_con),
       x = "AC1 (autocorrelation @ 1 s lag)", y = "Cells",
       fill = NULL) +
  theme_project() +
  theme(legend.position = "bottom")
print(p_ac1)

# --- Plot 3: n_beats distribution with the max_beats cap --------------------
p_nbeats <- contr_cell_params_qc %>%
  filter(!is.na(n_beats)) %>%
  ggplot(aes(x = n_beats, fill = signal_class)) +
  geom_histogram(binwidth = 1, colour = "white", alpha = 0.9) +
  geom_vline(xintercept = max_beats_con + 0.5,
             linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title = sprintf(
         "Contraction n_beats (cap = %d, cells with > cap are irregular)",
         max_beats_con),
       x = "n_beats", y = "Cells", fill = NULL) +
  theme_project() +
  theme(legend.position = "bottom")
print(p_nbeats)

# --- Plot 4: AC1 vs n_beats scatter, coloured by signal_class ---------------
p_scatter <- contr_cell_params_qc %>%
  filter(is.finite(ac1), !is.na(n_beats)) %>%
  ggplot(aes(x = ac1, y = n_beats, colour = signal_class)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = ac1_otsu_con,
             linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = max_beats_con + 0.5,
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller) +
  scale_colour_manual(values = class_colors) +
  labs(title = "Contraction: AC1 vs n_beats (QC decision plane)",
       x = "AC1", y = "n_beats", colour = NULL) +
  theme_project() +
  theme(legend.position = "bottom")
print(p_scatter)

# --- Plot 5: Example traces (good_paced and irregular) ---------------------
cat("  Plotting example traces...\n")
set.seed(42)

plot_example_traces <- function(indices, title_prefix) {
  if (length(indices) == 0) return(NULL)
  traces <- lapply(indices, function(i) {
    res <- contr_results[[i]]
    if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
    row <- contr_cell_params_qc[i, ]
    lbl <- paste0(row$cell_uid, " | ", row$Condition, " | ",
                  row$Cell_Line, " | ", row$signal_class,
                  "\nn=", row$n_beats,
                  "  AC1=", round(row$ac1, 3))
    tibble(
      label     = lbl,
      time_ms   = res$trace$time_ms,
      signal_bc = res$trace$signal_bc
    )
  })
  traces_df <- bind_rows(traces)
  if (nrow(traces_df) == 0) return(NULL)

  ggplot(traces_df, aes(x = time_ms, y = signal_bc)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ label, scales = "free", ncol = 3) +
    labs(title = title_prefix, x = "Time (ms)",
         y = "Contraction (raw stretch)") +
    theme_project() +
    theme(strip.text = element_text(size = 6))
}

gp_idx <- which(contr_cell_params_qc$signal_class == "good_paced")
if (length(gp_idx) > 6) gp_idx <- sample(gp_idx, 6)
p_pass <- plot_example_traces(gp_idx, "Example good_paced traces")
if (!is.null(p_pass)) print(p_pass)

irr_idx <- which(contr_cell_params_qc$signal_class == "irregular")
if (length(irr_idx) > 6) irr_idx <- sample(irr_idx, 6)
p_fail <- plot_example_traces(irr_idx, "Example irregular traces")
if (!is.null(p_fail)) print(p_fail)

## Close PDF ----------------------------------------------------------------
dev.off()
cat("  QC plots saved to:", pdf_path, "\n")

## Well-level QC summary table ---------------------------------------------
cat("  Generating well-level QC summary...\n")

# Reports the real QC outcome: fraction of cells classified as "good_paced"
# (AC1 >= Otsu + n_beats <= cap). `qc_pass` alone only marks whether the cell
# had detectable beats and is near-100% by construction, so it is not used
# here for the pass rate.
contr_well_qc <- contr_cell_params_qc %>%
  group_by(Plate, Well, well_uid, Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = sum(signal_class == "irregular"),
    pct_pass     = round(100 * n_good_paced / n_total, 1),
    .groups  = "drop"
  ) %>%
  arrange(pct_pass)

cat("  Well-level QC summary:\n")
print(as.data.frame(contr_well_qc), row.names = FALSE)

well_qc_path <- file.path(qc_dir, "contraction_well_qc_summary.csv")
write.csv(contr_well_qc, well_qc_path, row.names = FALSE)
cat("  Well QC summary saved to:", well_qc_path, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5b. QC Calibration gallery (stratified ~50 cells for labeling)
# ─────────────────────────────────────────────────────────────────────────────

cat("  Generating QC calibration gallery...\n")
set.seed(42)

# Stratify around the AC1 Otsu cutoff so the gallery focuses on the
# decision boundary (cells with AC1 just above / just below the threshold).
gp_cells  <- contr_cell_params_qc %>% filter(signal_class == "good_paced")
irr_cells <- contr_cell_params_qc %>% filter(signal_class == "irregular")

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
         !is.na(n_beats) & n_beats > max_beats_con) %>%
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

# Save companion CSV
calib_csv <- calib_cells %>%
  select(id, cell_uid, well_uid, Condition, Cell_Line,
         signal_class, stratum, ac1, n_beats, median_ibi_ms,
         snr_amp, rise_time_ms, decay_time_ms) %>%
  mutate(user_label = "")
calib_csv_path <- file.path(qc_dir,
  paste0("calibration_contraction_", analysis_date, ".csv"))
write.csv(calib_csv, calib_csv_path, row.names = FALSE)
cat("  Calibration CSV saved to:", calib_csv_path, "\n")

# Build trace data with peak markers
calib_trace_list <- lapply(seq_len(n_calib), function(j) {
  uid     <- calib_cells$cell_uid[j]
  row_idx <- which(singlecell$cell_uid == uid)[1]
  if (is.na(row_idx)) return(NULL)
  res <- contr_results[[row_idx]]
  if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
  ci <- calib_cells[j, ]
  label <- sprintf(
    "#%d | %s | %s | %s\nn=%s  AC1=%.3f  IBI=%.0f  snr=%.1f",
    ci$id, ci$cell_uid, ci$signal_class, ci$stratum,
    ifelse(is.na(ci$n_beats), "NA", as.character(ci$n_beats)),
    ifelse(is.na(ci$ac1), NA, ci$ac1),
    ifelse(is.na(ci$median_ibi_ms), NA, ci$median_ibi_ms),
    ifelse(is.na(ci$snr_amp), NA, ci$snr_amp)
  )
  # Trace data
  trace_df <- tibble(
    id           = ci$id,
    facet_label  = label,
    time_ms      = res$trace$time_ms,
    signal_bc    = res$trace$signal_bc,
    is_peak      = FALSE
  )
  # Add peak markers
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
  # Sort facets by ID number
  calib_trace_df$facet_label <- factor(
    calib_trace_df$facet_label,
    levels = unique(calib_trace_df$facet_label[
      order(calib_trace_df$id)]))

  traces_only <- calib_trace_df %>% filter(!is_peak)
  peaks_only  <- calib_trace_df %>% filter(is_peak)

  # Multi-page PDF: 6 traces per page (2 rows x 3 cols)
  calib_pdf_path <- file.path(qc_dir,
    paste0("calibration_contraction_", analysis_date, ".pdf"))
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
           y = "Contraction (raw stretch)") +
      theme_project() +
      theme(strip.text = element_text(size = 7))
    print(p)
  }
  dev.off()
  cat("  Calibration PDF saved to:", calib_pdf_path, "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 5c. Publication trace gallery
# ─────────────────────────────────────────────────────────────────────────────

cat("  Generating trace gallery...\n")

dir.create(file.path(main_dir, "Figures"), showWarnings = FALSE, recursive = TRUE)

## Pick 5 representative wells (top 3 per cell line, then trim to 5) ----------
gallery_well_pick <- contr_cell_params_qc %>%
  filter(signal_class == "good_paced") %>%
  count(well_uid, Cell_Line, name = "n_good_paced") %>%
  group_by(Cell_Line) %>%
  slice_max(n_good_paced, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  slice_max(n_good_paced, n = 5, with_ties = FALSE)

set.seed(123)

gallery_colors <- c("good_paced" = "#2166AC",
                    "irregular"  = "#D95F02")

## One PDF per well, 16 traces each (4 cols x 4 rows) ----------------------
for (w in seq_len(nrow(gallery_well_pick))) {
  wuid <- gallery_well_pick$well_uid[w]
  pool <- contr_cell_params_qc %>% filter(well_uid == wuid)

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
    res <- contr_results[[row_idx]]
    if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
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
    labs(title    = paste("Contraction gallery:", wuid),
         x = "Time (ms)", y = "Contraction (raw stretch)",
         color = "Signal class") +
    theme_project() +
    theme(strip.text      = element_text(size = 6),
          legend.position = "bottom")

  safe_wuid <- gsub("[^A-Za-z0-9_.]", "_", wuid)
  pdf_path  <- file.path(
    main_dir, "Figures",
    paste0("trace_gallery_contraction_", safe_wuid,
           "_", analysis_date, ".pdf")
  )
  cairo_pdf(pdf_path, width = 14, height = 12)
  print(p)
  dev.off()
  cat("  Gallery saved:", pdf_path, "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. Save outputs
# ─────────────────────────────────────────────────────────────────────────────

cat("Saving processed contraction data...\n")

## RDS with full cell params (including QC columns, before filtering)
rds_path <- file.path(processed_dir, "contraction_all_cells.rds")
saveRDS(contr_cell_params_qc, rds_path)
cat("  Saved:", rds_path, "\n")

## RDS with raw trace results (needed for cluster trace plots when skipping step 1)
results_rds_path <- file.path(processed_dir, "contraction_results_all.rds")
saveRDS(contr_results, results_rds_path)
cat("  Saved:", results_rds_path, "\n")

## Export to global environment for master script access
contr_cell_params_all  <<- contr_cell_params_qc
contr_results_all      <<- contr_results
contr_well_qc_summary  <<- contr_well_qc

cat("=== 1_process_contraction.R complete ===\n")
