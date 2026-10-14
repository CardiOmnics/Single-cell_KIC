# =============================================================================
# qc_plots_calcium.R
# Regenerate calcium QC plots from saved RDS files (use when run_calcium = FALSE).
#
# Expects in global env:
#   ca_cell_params_all, ca_results_all, main_dir, analysis_date,
#   theme_project(), otsu_threshold()
# =============================================================================

cat("=== qc_plots_calcium.R ===\n")

# ---- Setup ------------------------------------------------------------------
qc_dir        <- file.path(main_dir, "QC")
processed_dir <- file.path(main_dir, "processed")
dir.create(qc_dir,        showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Aliases matching the names used in the original processing script
ca_cell_params_qc <- ca_cell_params_all
ca_results        <- ca_results_all

# Recompute thresholds from the stored ac1 column
max_beats_ca <- 10
ac1_otsu_ca  <- otsu_threshold(ca_cell_params_qc$ac1)
cat(sprintf("  AC1 Otsu threshold (recomputed): %.4f\n", ac1_otsu_ca))

# Always recompute the well QC gate from current data so this script is
# safe to run standalone (avoids stale ca_well_qc_gate from a previous run).
if (!exists("min_good_paced")) min_good_paced <- 10
ca_well_qc_gate <- ca_cell_params_qc %>%
  group_by(well_uid, Plate, Well, Condition, Cell_Line) %>%
  summarise(
    n_total        = n(),
    n_good_paced   = sum(signal_class == "good_paced"),
    n_irregular    = sum(signal_class == "irregular"),
    n_manually_excluded = sum(!qc_pass),
    pct_good_paced = round(100 * n_good_paced / n_total, 1),
    pct_irregular  = round(100 * n_irregular  / n_total, 1),
    .groups        = "drop"
  ) %>%
  mutate(
    well_qc_pass   = n_good_paced >= min_good_paced,
    qc_fail_reason = case_when(
      !well_qc_pass ~ sprintf("n_good_paced = %d (< %d)",
                               n_good_paced, min_good_paced),
      n_manually_excluded == n_total ~ "manually excluded",
      TRUE ~ NA_character_
    ),
    modality       = "Calcium"
  )
cat(sprintf("  Well QC gate (recomputed): %d wells, %d failing\n",
            nrow(ca_well_qc_gate), sum(!ca_well_qc_gate$well_qc_pass)))

# ---- Main QC plots PDF ------------------------------------------------------
pdf_path <- file.path(qc_dir, "calcium_qc_plots.pdf")
cairo_pdf(pdf_path, width = 16, height = 10, onefile = TRUE)

# Per-well trace overlay
cat("    Per-well trace overlay...\n")
trace_list <- lapply(seq_along(ca_results), function(i) {
  res <- ca_results[[i]]
  if (is.null(res) || is.null(res$trace)) return(NULL)
  res$trace %>%
    mutate(
      cell_uid  = ca_cell_params_qc$cell_uid[i],
      well_uid  = ca_cell_params_qc$well_uid[i],
      Condition = ca_cell_params_qc$Condition[i],
      Cell_Line = ca_cell_params_qc$Cell_Line[i]
    )
})
trace_df <- bind_rows(trace_list)

if (nrow(trace_df) > 0) {
  wells <- sort(unique(trace_df$well_uid))
  n_per_page <- 20
  n_pages <- ceiling(length(wells) / n_per_page)
  for (pg in seq_len(n_pages)) {
    wells_page <- wells[seq((pg - 1) * n_per_page + 1, min(pg * n_per_page, length(wells)))]
    p_traces <- trace_df %>%
      filter(well_uid %in% wells_page) %>%
      ggplot(aes(x = time_ms, y = signal_bc, group = cell_uid)) +
      geom_line(alpha = 0.15, linewidth = 0.2) +
      facet_wrap(~ well_uid, scales = "free_y") +
      labs(title = sprintf("Calcium traces per well (page %d/%d)", pg, n_pages),
           x = "Time (ms)", y = "Calcium signal (CyteSeer Trimmed)") +
      theme_project()
    print(p_traces)
  }
}

# AC1 distribution
cat("    Plotting AC1 distribution...\n")
class_colors <- c(good_paced = "#0072B2", irregular = "#D55E00")

p_ac1 <- ca_cell_params_qc %>%
  filter(is.finite(ac1)) %>%
  ggplot(aes(x = ac1, fill = signal_class)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
  geom_vline(xintercept = ac1_otsu_ca, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title   = sprintf("Calcium AC1 distribution (Otsu threshold = %.3f)", ac1_otsu_ca),
       subtitle = sprintf("good_paced also requires n_beats <= %d", max_beats_ca),
       x = "AC1 (autocorrelation @ 1 s lag)", y = "Cells", fill = NULL) +
  theme_project() + theme(legend.position = "bottom")
print(p_ac1)

# n_beats distribution
p_nbeats <- ca_cell_params_qc %>%
  filter(!is.na(n_beats)) %>%
  ggplot(aes(x = n_beats, fill = signal_class)) +
  geom_histogram(binwidth = 1, colour = "white", alpha = 0.9) +
  geom_vline(xintercept = max_beats_ca + 0.5, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title = sprintf("Calcium n_beats (cap = %d)", max_beats_ca),
       x = "n_beats", y = "Cells", fill = NULL) +
  theme_project() + theme(legend.position = "bottom")
print(p_nbeats)

# AC1 vs n_beats scatter
p_scatter <- ca_cell_params_qc %>%
  filter(is.finite(ac1), !is.na(n_beats)) %>%
  ggplot(aes(x = ac1, y = n_beats, colour = signal_class)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = ac1_otsu_ca, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = max_beats_ca + 0.5, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller) +
  scale_colour_manual(values = class_colors) +
  labs(title = "Calcium: AC1 vs n_beats (QC decision plane)",
       x = "AC1", y = "n_beats", colour = NULL) +
  theme_project() + theme(legend.position = "bottom")
print(p_scatter)

# Morphology diagnostic -------------------------------------------------------
# Independent of the Ca signal: is the segmented-cell morphology unimodal
# (consistent with a homogeneous CM-dominated population) or bimodal
# (suggesting a second cell type)? Also: does low-AC1 cluster at a distinct
# morphology, or span the same range as high-AC1? If it spans, low AC1 is
# mostly CMs with imaging issues, not non-CMs.
cat("    Plotting morphology diagnostic...\n")

morph_df <- ca_cell_params_qc %>%
  filter(is.finite(Cell.Area), is.finite(Nucleus.Area),
         Cell.Area > 0, Nucleus.Area > 0) %>%
  mutate(ac1_bin = case_when(
    !is.finite(ac1)            ~ "NA",
    ac1 >= ac1_otsu_ca         ~ "AC1 >= Otsu",
    TRUE                       ~ "AC1 <  Otsu"
  ))

if (nrow(morph_df) > 0) {
  # Cell area distribution per cell line — check for bimodality
  p_area <- morph_df %>%
    ggplot(aes(x = Cell.Area)) +
    geom_density(fill = "grey60", colour = NA, alpha = 0.7) +
    facet_wrap(~ Cell_Line, scales = "free_y",
               labeller = cell_line_fig_labeller) +
    scale_x_log10() +
    labs(title = "Cell area distribution (independent of Ca signal)",
         subtitle = "Unimodal = consistent with homogeneous CM population",
         x = expression("Cell area ("*mu*m^2*", log scale)"),
         y = "Density") +
    theme_project()
  print(p_area)

  # NC ratio distribution per cell line
  p_nc <- morph_df %>%
    ggplot(aes(x = NC_ratio)) +
    geom_density(fill = "grey60", colour = NA, alpha = 0.7) +
    facet_wrap(~ Cell_Line, scales = "free_y",
               labeller = cell_line_fig_labeller) +
    labs(title = "Nucleus-to-cell area ratio",
         subtitle = "Non-CMs (e.g. fibroblasts) tend to have higher NC ratio",
         x = "Nucleus area / Cell area",
         y = "Density") +
    theme_project()
  print(p_nc)

  # AC1-bin vs morphology — does low AC1 cluster at distinct morphology?
  ac1_morph_colors <- c("AC1 >= Otsu" = "#0072B2",
                        "AC1 <  Otsu" = "#D55E00",
                        "NA"          = "grey70")
  p_morph_ac1 <- morph_df %>%
    ggplot(aes(x = Cell.Area, y = NC_ratio, colour = ac1_bin)) +
    geom_point(alpha = 0.4, size = 0.6) +
    facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller) +
    scale_x_log10() +
    scale_colour_manual(values = ac1_morph_colors) +
    labs(title = "Morphology vs AC1 status",
         subtitle = paste(
           "If low-AC1 cells span the high-AC1 morphology range,",
           "low AC1 is mostly CMs with imaging issues, not non-CMs"),
         x = expression("Cell area ("*mu*m^2*", log scale)"),
         y = "NC ratio", colour = NULL) +
    theme_project() +
    theme(legend.position = "bottom")
  print(p_morph_ac1)
}

# Example traces
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
    res$trace %>% mutate(facet_label = label, example_group = cls)
  }))
}
example_trace_df <- bind_rows(sample_examples("good_paced"), sample_examples("irregular"))
if (nrow(example_trace_df) > 0) {
  for (grp in unique(example_trace_df$example_group)) {
    grp_df <- example_trace_df %>% filter(example_group == grp)
    if (nrow(grp_df) == 0) next
    p_ex <- ggplot(grp_df, aes(x = time_ms, y = signal_bc)) +
      geom_line(linewidth = 0.3) +
      facet_wrap(~ facet_label, scales = "free_y") +
      labs(title = sprintf("Calcium: example traces — %s", grp),
           x = "Time (ms)", y = "Signal") +
      theme_project() + theme(strip.text = element_text(size = 7))
    print(p_ex)
  }
}
dev.off()
cat(sprintf("  Saved QC plots: %s\n", pdf_path))

# ---- Well-level QC summary --------------------------------------------------
ca_well_qc <- ca_cell_params_qc %>%
  group_by(Plate, Well, well_uid, Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = sum(signal_class == "irregular"),
    pct_pass     = round(100 * n_good_paced / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(pct_pass)
well_qc_path <- file.path(qc_dir, "calcium_well_qc_summary.csv")
write.csv(ca_well_qc, well_qc_path, row.names = FALSE)
cat(sprintf("  Saved well QC summary: %s\n", well_qc_path))

# ---- QC calibration gallery -------------------------------------------------
cat("  Generating QC calibration gallery...\n")
set.seed(42)
gp_cells  <- ca_cell_params_qc %>% filter(signal_class == "good_paced")
irr_cells <- ca_cell_params_qc %>% filter(signal_class == "irregular")
calib_cells <- bind_rows(
  gp_cells  %>% filter(is.finite(ac1)) %>% slice_min(ac1, n = 15, with_ties = FALSE) %>% mutate(stratum = "gp_just_above_otsu"),
  irr_cells %>% filter(is.finite(ac1)) %>% slice_max(ac1, n = 10, with_ties = FALSE) %>% mutate(stratum = "irr_just_below_otsu"),
  irr_cells %>% filter(is.finite(ac1), !is.na(n_beats) & n_beats > max_beats_ca) %>% slice_sample(n = min(10, nrow(.))) %>% mutate(stratum = "irr_nbeats_cap"),
  gp_cells  %>% filter(!cell_uid %in% (gp_cells  %>% filter(is.finite(ac1)) %>% slice_min(ac1, n = 15, with_ties = FALSE))$cell_uid) %>% slice_sample(n = min(10, nrow(.))) %>% mutate(stratum = "gp_random"),
  irr_cells %>% slice_sample(n = min(5, nrow(.))) %>% mutate(stratum = "irr_random")
) %>% distinct(cell_uid, .keep_all = TRUE) %>% mutate(id = row_number())
n_calib <- nrow(calib_cells)
cat("  Calibration sample:", n_calib, "cells\n")

calib_csv <- calib_cells %>%
  select(id, cell_uid, well_uid, Condition, Cell_Line,
         signal_class, stratum, ac1, n_beats, median_ibi_ms,
         snr_amp, rise_time_ms, decay_time_ms,
         mean_amp, CTD50_ms, CTD90_ms) %>%
  mutate(user_label = "")
write.csv(calib_csv,
          file.path(qc_dir, paste0("calibration_calcium_", analysis_date, ".csv")),
          row.names = FALSE)

calib_trace_list <- lapply(seq_len(n_calib), function(j) {
  uid     <- calib_cells$cell_uid[j]
  row_idx <- which(ca_cell_params_qc$cell_uid == uid)[1]
  if (is.na(row_idx)) return(NULL)
  res <- ca_results[[row_idx]]
  if (is.null(res) || is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
  ci <- calib_cells[j, ]
  label <- sprintf("#%d | %s | %s | %s\nn=%s  AC1=%.3f  IBI=%.0f  snr=%.1f",
    ci$id, ci$cell_uid, ci$signal_class, ci$stratum,
    ifelse(is.na(ci$n_beats), "NA", as.character(ci$n_beats)),
    ifelse(is.na(ci$ac1), NA, ci$ac1),
    ifelse(is.na(ci$median_ibi_ms), NA, ci$median_ibi_ms),
    ifelse(is.na(ci$snr_amp), NA, ci$snr_amp))
  trace_df <- tibble(id = ci$id, facet_label = label,
                     time_ms = res$trace$time_ms, signal_bc = res$trace$signal_bc, is_peak = FALSE)
  if (!is.null(res$peaks) && nrow(res$peaks) > 0)
    trace_df <- bind_rows(trace_df,
      tibble(id = ci$id, facet_label = label,
             time_ms = res$peaks$peak_time, signal_bc = res$peaks$peak_amp, is_peak = TRUE))
  trace_df
})
calib_trace_df <- bind_rows(calib_trace_list)
if (nrow(calib_trace_df) > 0) {
  calib_trace_df$facet_label <- factor(calib_trace_df$facet_label,
    levels = unique(calib_trace_df$facet_label[order(calib_trace_df$id)]))
  calib_pdf_path <- file.path(qc_dir, paste0("calibration_calcium_", analysis_date, ".pdf"))
  page_ids <- split(sort(unique(calib_cells$id)),
                    ceiling(seq_along(unique(calib_cells$id)) / 6))
  cairo_pdf(calib_pdf_path, width = 14, height = 10, onefile = TRUE)
  for (pg in page_ids) {
    pg_labels <- levels(calib_trace_df$facet_label)[pg]
    p <- ggplot(calib_trace_df %>% filter(!is_peak, facet_label %in% pg_labels),
                aes(x = time_ms, y = signal_bc)) +
      geom_line(linewidth = 0.4, color = "grey30") +
      geom_point(data = calib_trace_df %>% filter(is_peak, facet_label %in% pg_labels),
                 color = "red", size = 2) +
      facet_wrap(~ facet_label, scales = "free", ncol = 3, nrow = 2) +
      labs(x = "Time (ms)", y = "Calcium signal (CyteSeer Trimmed)") +
      theme_project() + theme(strip.text = element_text(size = 7))
    print(p)
  }
  dev.off()
  cat("  Calibration PDF saved to:", calib_pdf_path, "\n")
}

# ---- Publication trace gallery ----------------------------------------------
cat("  Generating trace gallery...\n")
if (!exists("expl_dir"))
  expl_dir <- file.path(main_dir, "Figures", "Exploratory")
gallery_dir <- file.path(expl_dir, "trace_galleries")
dir.create(gallery_dir, showWarnings = FALSE, recursive = TRUE)
gallery_colors <- c("good_paced" = "#2166AC", "irregular" = "#D95F02")
gallery_well_pick <- ca_cell_params_qc %>%
  filter(signal_class == "good_paced") %>%
  count(well_uid, Cell_Line, name = "n_good_paced") %>%
  group_by(Cell_Line) %>% slice_max(n_good_paced, n = 3, with_ties = FALSE) %>%
  ungroup() %>% slice_max(n_good_paced, n = 5, with_ties = FALSE)
set.seed(123)
for (w in seq_len(nrow(gallery_well_pick))) {
  wuid <- gallery_well_pick$well_uid[w]
  pool <- ca_cell_params_qc %>% filter(well_uid == wuid)
  picked <- bind_rows(
    if (nrow(pool %>% filter(signal_class == "good_paced")) > 0)
      slice_sample(pool %>% filter(signal_class == "good_paced"), n = min(12, nrow(pool %>% filter(signal_class == "good_paced")))) else NULL,
    if (nrow(pool %>% filter(signal_class == "irregular")) > 0)
      slice_sample(pool %>% filter(signal_class == "irregular"), n = min(4, nrow(pool %>% filter(signal_class == "irregular")))) else NULL
  )
  if (nrow(picked) > 16) picked <- slice_sample(picked, n = 16)
  trace_list <- lapply(seq_len(nrow(picked)), function(j) {
    uid     <- picked$cell_uid[j]
    row_idx <- which(ca_cell_params_qc$cell_uid == uid)[1]
    if (is.na(row_idx)) return(NULL)
    res <- ca_results[[row_idx]]
    if (is.null(res) || is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
    ci <- picked[j, ]
    label <- sprintf("%s\n%s | n=%s  AC1=%.3f  medIBI=%.0f  snr=%.1f",
      uid, ci$signal_class,
      ifelse(is.na(ci$n_beats), "NA", as.character(ci$n_beats)),
      ifelse(is.na(ci$ac1), NA, ci$ac1),
      ifelse(is.na(ci$median_ibi_ms), NA, ci$median_ibi_ms),
      ifelse(is.na(ci$snr_amp), NA, ci$snr_amp))
    tibble(facet_label = label, signal_class = ci$signal_class,
           time_ms = res$trace$time_ms, signal_bc = res$trace$signal_bc)
  })
  trace_df <- bind_rows(trace_list)
  if (nrow(trace_df) == 0) next
  p <- ggplot(trace_df, aes(x = time_ms, y = signal_bc, color = signal_class)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ facet_label, scales = "free", ncol = 4) +
    scale_color_manual(values = gallery_colors) +
    labs(title = paste("Calcium gallery:", wuid),
         x = "Time (ms)", y = "Calcium signal (CyteSeer Trimmed)", color = "Signal class") +
    theme_project() +
    theme(strip.text = element_text(size = 6), legend.position = "bottom")
  safe_wuid <- gsub("[^A-Za-z0-9_.]", "_", wuid)
  pdf_path  <- file.path(gallery_dir,
    paste0("trace_gallery_calcium_", safe_wuid, "_", analysis_date, ".pdf"))
  cairo_pdf(pdf_path, width = 14, height = 12)
  print(p)
  dev.off()
  cat("  Gallery saved:", pdf_path, "\n")
}

# ---- Well QC gallery: trace overlay per well, colored by signal class ----------
# Irregular traces = red, good_paced = black. Failed well strip text = red.
cat("  Generating well QC gallery...\n")

  # Rebuild trace data from ca_results (trace_df may have been overwritten
# by the publication gallery loop above)
trace_all <- bind_rows(lapply(seq_along(ca_results), function(i) {
  res <- ca_results[[i]]
  if (is.null(res) || is.null(res$trace)) return(NULL)
  res$trace %>%
    mutate(
      cell_uid  = ca_cell_params_qc$cell_uid[i],
      well_uid  = ca_cell_params_qc$well_uid[i]
    )
}))

if (exists("ca_well_qc_gate") && nrow(trace_all) > 0) {
  # Add signal_class to trace data
  cell_classes <- ca_cell_params_qc %>%
    select(cell_uid, signal_class)
  trace_qc <- trace_all %>%
    left_join(cell_classes, by = "cell_uid")

  # Build facet labels with QC status — embed HTML color for per-panel styling
  if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext")
  library(ggtext)

  well_labels <- ca_well_qc_gate %>%
    mutate(
      facet_label = ifelse(
        well_qc_pass,
        sprintf("%s | %s | %s<br>n_gp=%d (%s%%)",
                well_uid, Condition, Cell_Line,
                n_good_paced, pct_good_paced),
        sprintf("<span style='color:red'>%s | %s | %s<br>FAILED: %s</span>",
                well_uid, Condition, Cell_Line,
                qc_fail_reason)
      )
    ) %>%
    select(well_uid, facet_label, well_qc_pass)

  trace_qc <- trace_qc %>%
    left_join(well_labels, by = "well_uid") %>%
    filter(!is.na(facet_label))

  # Sort: failed wells first, then by well_uid
  label_order <- well_labels %>%
    arrange(well_qc_pass, well_uid) %>%
    pull(facet_label)
  trace_qc$facet_label <- factor(trace_qc$facet_label, levels = label_order)

  wells_per_page <- 12
  all_labels <- levels(trace_qc$facet_label)
  n_pages <- ceiling(length(all_labels) / wells_per_page)

  well_qc_pdf <- file.path(qc_dir,
    paste0("well_qc_calcium_", analysis_date, ".pdf"))
  cairo_pdf(well_qc_pdf, width = 16, height = 10, onefile = TRUE)

  for (pg in seq_len(n_pages)) {
    idx_start <- (pg - 1) * wells_per_page + 1
    idx_end   <- min(pg * wells_per_page, length(all_labels))
    page_labels <- all_labels[idx_start:idx_end]

    page_traces <- trace_qc %>% filter(facet_label %in% page_labels)

    p <- ggplot(page_traces,
                aes(x = time_ms, y = signal_bc,
                    group = cell_uid, colour = signal_class)) +
      geom_line(alpha = 0.3, linewidth = 0.2) +
      facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
      scale_colour_manual(
        values = c(good_paced = "black", irregular = "red"),
        labels = c(good_paced = "Good paced", irregular = "Irregular"),
        name = NULL) +
      labs(title = sprintf("Calcium: well QC overview (page %d/%d)", pg, n_pages),
           subtitle = sprintf("Red title = failed well QC (< %d good_paced cells) | Red traces = irregular",
                              min_good_paced),
           x = "Time (ms)", y = "Calcium signal (CyteSeer Trimmed)") +
      theme_project() +
      theme(strip.text = ggtext::element_markdown(size = 6),
            legend.position = "bottom")
    print(p)
  }

  dev.off()
  cat("  Well QC gallery saved to:", well_qc_pdf, "\n")
} else {
  cat("  Skipped well QC gallery (no ca_well_qc_gate or traces)\n")
}

cat("=== qc_plots_calcium.R complete ===\n")
