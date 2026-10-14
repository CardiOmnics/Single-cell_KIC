# =============================================================================
# qc_plots_contraction.R
# Regenerate contraction QC plots from saved RDS files (use when run_contraction = FALSE).
#
# Expects in global env:
#   contr_cell_params_all, contr_results_all, main_dir, analysis_date,
#   theme_project(), otsu_threshold()
# =============================================================================

cat("=== qc_plots_contraction.R ===\n")

# ---- Setup ------------------------------------------------------------------
qc_dir        <- file.path(main_dir, "QC")
processed_dir <- file.path(main_dir, "processed")
dir.create(qc_dir,        showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Aliases matching the names used in the original processing script
contr_cell_params_qc <- contr_cell_params_all
contr_results        <- contr_results_all

# Recompute thresholds from the stored ac1 column
max_beats_con <- 5
ac1_otsu_con  <- otsu_threshold(contr_cell_params_qc$ac1)
cat(sprintf("  AC1 Otsu threshold (recomputed): %.4f\n", ac1_otsu_con))

# Always recompute the well QC gate from current data so this script is
# safe to run standalone (avoids stale contr_well_qc_gate from a previous run).
if (!exists("min_good_paced")) min_good_paced <- 10
contr_well_qc_gate <- contr_cell_params_qc %>%
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
    modality       = "Contraction"
  )
cat(sprintf("  Well QC gate (recomputed): %d wells, %d failing\n",
            nrow(contr_well_qc_gate), sum(!contr_well_qc_gate$well_qc_pass)))

# ---- Main QC plots PDF ------------------------------------------------------
pdf_path <- file.path(qc_dir, "contraction_qc_plots.pdf")
cairo_pdf(pdf_path, width = 14, height = 10, onefile = TRUE)

# Per-well trace overlay
cat("  Building per-well trace overlay...\n")
trace_list <- lapply(seq_along(contr_results), function(i) {
  res <- contr_results[[i]]
  if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
  tibble(
    cell_uid  = contr_cell_params_qc$cell_uid[i],
    well_uid  = contr_cell_params_qc$well_uid[i],
    Condition = as.character(contr_cell_params_qc$Condition[i]),
    Cell_Line = as.character(contr_cell_params_qc$Cell_Line[i]),
    time_ms   = res$trace$time_ms,
    signal_bc = res$trace$signal_bc
  )
})
trace_long <- bind_rows(trace_list)

if (nrow(trace_long) > 0) {
  wells_all      <- sort(unique(trace_long$well_uid))
  wells_per_page <- 12
  n_pages        <- ceiling(length(wells_all) / wells_per_page)
  for (pg in seq_len(n_pages)) {
    wells_sub <- wells_all[seq((pg - 1) * wells_per_page + 1, min(pg * wells_per_page, length(wells_all)))]
    p_traces <- trace_long %>%
      filter(well_uid %in% wells_sub) %>%
      ggplot(aes(x = time_ms, y = signal_bc, group = cell_uid)) +
      geom_line(alpha = 0.15, linewidth = 0.2) +
      facet_wrap(~ well_uid, scales = "free_y", ncol = 4) +
      labs(title = paste0("Contraction trace overlay (page ", pg, "/", n_pages, ")"),
           x = "Time (ms)", y = "Contraction (raw stretch)") +
      theme_project() + theme(strip.text = element_text(size = 7))
    print(p_traces)
  }
} else {
  plot.new(); text(0.5, 0.5, "No contraction traces available", cex = 1.2)
}

# AC1 distribution
cat("  Plotting AC1 distribution...\n")
class_colors <- c(good_paced = "#2166AC", irregular = "#D95F02")
p_ac1 <- contr_cell_params_qc %>%
  filter(is.finite(ac1)) %>%
  ggplot(aes(x = ac1, fill = signal_class)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
  geom_vline(xintercept = ac1_otsu_con, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title   = sprintf("Contraction AC1 distribution (Otsu threshold = %.3f)", ac1_otsu_con),
       subtitle = sprintf("good_paced also requires n_beats <= %d", max_beats_con),
       x = "AC1 (autocorrelation @ 1 s lag)", y = "Cells", fill = NULL) +
  theme_project() + theme(legend.position = "bottom")
print(p_ac1)

# n_beats distribution
p_nbeats <- contr_cell_params_qc %>%
  filter(!is.na(n_beats)) %>%
  ggplot(aes(x = n_beats, fill = signal_class)) +
  geom_histogram(binwidth = 1, colour = "white", alpha = 0.9) +
  geom_vline(xintercept = max_beats_con + 0.5, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ Cell_Line, scales = "free_y", labeller = cell_line_fig_labeller) +
  scale_fill_manual(values = class_colors) +
  labs(title = sprintf("Contraction n_beats (cap = %d)", max_beats_con),
       x = "n_beats", y = "Cells", fill = NULL) +
  theme_project() + theme(legend.position = "bottom")
print(p_nbeats)

# AC1 vs n_beats scatter
p_scatter <- contr_cell_params_qc %>%
  filter(is.finite(ac1), !is.na(n_beats)) %>%
  ggplot(aes(x = ac1, y = n_beats, colour = signal_class)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = ac1_otsu_con, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = max_beats_con + 0.5, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller) +
  scale_colour_manual(values = class_colors) +
  labs(title = "Contraction: AC1 vs n_beats (QC decision plane)",
       x = "AC1", y = "n_beats", colour = NULL) +
  theme_project() + theme(legend.position = "bottom")
print(p_scatter)

# Example traces
cat("  Plotting example traces...\n")
set.seed(42)
plot_example_traces <- function(indices, title_prefix) {
  if (length(indices) == 0) return(NULL)
  traces <- lapply(indices, function(i) {
    res <- contr_results[[i]]
    if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
    row <- contr_cell_params_qc[i, ]
    lbl <- paste0(row$cell_uid, " | ", row$Condition, " | ", row$Cell_Line, " | ",
                  row$signal_class, "\nn=", row$n_beats, "  AC1=", round(row$ac1, 3))
    tibble(label = lbl, time_ms = res$trace$time_ms, signal_bc = res$trace$signal_bc)
  })
  traces_df <- bind_rows(traces)
  if (nrow(traces_df) == 0) return(NULL)
  ggplot(traces_df, aes(x = time_ms, y = signal_bc)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ label, scales = "free", ncol = 3) +
    labs(title = title_prefix, x = "Time (ms)", y = "Contraction (raw stretch)") +
    theme_project() + theme(strip.text = element_text(size = 6))
}
gp_idx  <- which(contr_cell_params_qc$signal_class == "good_paced")
if (length(gp_idx)  > 6) gp_idx  <- sample(gp_idx,  6)
irr_idx <- which(contr_cell_params_qc$signal_class == "irregular")
if (length(irr_idx) > 6) irr_idx <- sample(irr_idx, 6)
p_pass <- plot_example_traces(gp_idx,  "Example good_paced traces")
p_fail <- plot_example_traces(irr_idx, "Example irregular traces")
if (!is.null(p_pass)) print(p_pass)
if (!is.null(p_fail)) print(p_fail)

dev.off()
cat("  QC plots saved to:", pdf_path, "\n")

# ---- Well-level QC summary --------------------------------------------------
contr_well_qc <- contr_cell_params_qc %>%
  group_by(Plate, Well, well_uid, Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = sum(signal_class == "irregular"),
    pct_pass     = round(100 * n_good_paced / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(pct_pass)
cat("  Well-level QC summary:\n")
print(as.data.frame(contr_well_qc), row.names = FALSE)
well_qc_path <- file.path(qc_dir, "contraction_well_qc_summary.csv")
write.csv(contr_well_qc, well_qc_path, row.names = FALSE)
cat("  Well QC summary saved to:", well_qc_path, "\n")

# ---- QC calibration gallery -------------------------------------------------
cat("  Generating QC calibration gallery...\n")
set.seed(42)
gp_cells  <- contr_cell_params_qc %>% filter(signal_class == "good_paced")
irr_cells <- contr_cell_params_qc %>% filter(signal_class == "irregular")
calib_cells <- bind_rows(
  gp_cells  %>% filter(is.finite(ac1)) %>% slice_min(ac1, n = 15, with_ties = FALSE) %>% mutate(stratum = "gp_just_above_otsu"),
  irr_cells %>% filter(is.finite(ac1)) %>% slice_max(ac1, n = 10, with_ties = FALSE) %>% mutate(stratum = "irr_just_below_otsu"),
  irr_cells %>% filter(is.finite(ac1), !is.na(n_beats) & n_beats > max_beats_con) %>% slice_sample(n = min(10, nrow(.))) %>% mutate(stratum = "irr_nbeats_cap"),
  gp_cells  %>% filter(!cell_uid %in% (gp_cells %>% filter(is.finite(ac1)) %>% slice_min(ac1, n = 15, with_ties = FALSE))$cell_uid) %>% slice_sample(n = min(10, nrow(.))) %>% mutate(stratum = "gp_random"),
  irr_cells %>% slice_sample(n = min(5, nrow(.))) %>% mutate(stratum = "irr_random")
) %>% distinct(cell_uid, .keep_all = TRUE) %>% mutate(id = row_number())
n_calib <- nrow(calib_cells)
cat("  Calibration sample:", n_calib, "cells\n")

calib_csv <- calib_cells %>%
  select(id, cell_uid, well_uid, Condition, Cell_Line,
         signal_class, stratum, ac1, n_beats, median_ibi_ms,
         snr_amp, rise_time_ms, decay_time_ms) %>%
  mutate(user_label = "")
write.csv(calib_csv,
          file.path(qc_dir, paste0("calibration_contraction_", analysis_date, ".csv")),
          row.names = FALSE)

calib_trace_list <- lapply(seq_len(n_calib), function(j) {
  uid     <- calib_cells$cell_uid[j]
  row_idx <- which(contr_cell_params_qc$cell_uid == uid)[1]
  if (is.na(row_idx)) return(NULL)
  res <- contr_results[[row_idx]]
  if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
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
  calib_pdf_path <- file.path(qc_dir, paste0("calibration_contraction_", analysis_date, ".pdf"))
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
      labs(x = "Time (ms)", y = "Contraction (raw stretch)") +
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
gallery_well_pick <- contr_cell_params_qc %>%
  filter(signal_class == "good_paced") %>%
  count(well_uid, Cell_Line, name = "n_good_paced") %>%
  group_by(Cell_Line) %>% slice_max(n_good_paced, n = 3, with_ties = FALSE) %>%
  ungroup() %>% slice_max(n_good_paced, n = 5, with_ties = FALSE)
set.seed(123)
for (w in seq_len(nrow(gallery_well_pick))) {
  wuid <- gallery_well_pick$well_uid[w]
  pool <- contr_cell_params_qc %>% filter(well_uid == wuid)
  picked <- bind_rows(
    if (nrow(pool %>% filter(signal_class == "good_paced")) > 0)
      slice_sample(pool %>% filter(signal_class == "good_paced"), n = min(12, nrow(pool %>% filter(signal_class == "good_paced")))) else NULL,
    if (nrow(pool %>% filter(signal_class == "irregular")) > 0)
      slice_sample(pool %>% filter(signal_class == "irregular"), n = min(4, nrow(pool %>% filter(signal_class == "irregular")))) else NULL
  )
  if (nrow(picked) > 16) picked <- slice_sample(picked, n = 16)
  trace_list <- lapply(seq_len(nrow(picked)), function(j) {
    uid     <- picked$cell_uid[j]
    row_idx <- which(contr_cell_params_qc$cell_uid == uid)[1]
    if (is.na(row_idx)) return(NULL)
    res <- contr_results[[row_idx]]
    if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
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
    labs(title = paste("Contraction gallery:", wuid),
         x = "Time (ms)", y = "Contraction (raw stretch)", color = "Signal class") +
    theme_project() +
    theme(strip.text = element_text(size = 6), legend.position = "bottom")
  safe_wuid <- gsub("[^A-Za-z0-9_.]", "_", wuid)
  pdf_path  <- file.path(gallery_dir,
    paste0("trace_gallery_contraction_", safe_wuid, "_", analysis_date, ".pdf"))
  cairo_pdf(pdf_path, width = 14, height = 12)
  print(p)
  dev.off()
  cat("  Gallery saved:", pdf_path, "\n")
}

# ---- Well QC gallery: trace overlay per well, colored by signal class ----------
# Irregular traces = red, good_paced = black. Failed well strip text = red.
cat("  Generating well QC gallery...\n")

if (exists("contr_well_qc_gate") && nrow(trace_long) > 0) {
  # Add signal_class to trace data
  cell_classes <- contr_cell_params_qc %>%
    select(cell_uid, signal_class)
  trace_qc <- trace_long %>%
    left_join(cell_classes, by = "cell_uid")

  # Build facet labels with QC status — embed HTML color for per-panel styling
  if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext")
  library(ggtext)

  well_labels <- contr_well_qc_gate %>%
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
    paste0("well_qc_contraction_", analysis_date, ".pdf"))
  cairo_pdf(well_qc_pdf, width = 14, height = 10, onefile = TRUE)

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
      labs(title = sprintf("Contraction: well QC overview (page %d/%d)", pg, n_pages),
           subtitle = sprintf("Red title = failed well QC (< %d good_paced cells) | Red traces = irregular",
                              min_good_paced),
           x = "Time (ms)", y = "Contraction (raw stretch)") +
      theme_project() +
      theme(strip.text = ggtext::element_markdown(size = 6),
            legend.position = "bottom")
    print(p)
  }

  dev.off()
  cat("  Well QC gallery saved to:", well_qc_pdf, "\n")
} else {
  cat("  Skipped well QC gallery (no contr_well_qc_gate or traces)\n")
}

cat("=== qc_plots_contraction.R complete ===\n")
