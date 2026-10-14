# =============================================================================
# 7_figure1_trace_illustration.R
# Figure 1 (trace-extraction panels): 8 individual single-cell traces per
# modality (calcium / contraction), small-multiples layout (4 x 2 grid).
# No axis labels — each cell gets a labelled sub-panel. The rest of Figure 1
# is a hand-drawn schematic.
#
# The well is drawn at random among wells that passed well-level QC in both
# modalities and have at least n_cells good-paced cells in each, with a fixed
# seed so the same well is chosen on every run.
#
# Uses ca_cell_params_all / contr_cell_params_all (well QC applied) and
# ca_results_all / contr_results_all from run_pipeline.R.
# Output:   3_output/Figures/Publication/figure1_trace_illustration/
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("\n========== 7_figure1_trace_illustration.R ==========\n")

if (!exists("main_dir"))      main_dir      <- "../3_output"
if (!exists("analysis_date")) analysis_date <- format(Sys.Date(), "%Y-%m-%d")
if (!exists("pub_dir"))       pub_dir       <- file.path(main_dir, "Figures", "Publication")
n_cells       <- 8
n_col         <- 4          # columns in facet grid  → 4 x 2
trace_window  <- c(0, 5000) # ms

# ---- Data (well-level QC already applied to qc_pass by run_pipeline.R) ----
ca_params   <- ca_cell_params_all
ca_results  <- ca_results_all
con_params  <- contr_cell_params_all
con_results <- contr_results_all

# Wells with >= n_cells good-paced, QC-passing cells in both modalities
n_good_per_well <- function(params) {
  params %>%
    filter(qc_pass, signal_class == "good_paced") %>%
    count(well_uid, name = "n_good")
}
eligible <- inner_join(n_good_per_well(ca_params), n_good_per_well(con_params),
                       by = "well_uid", suffix = c("_ca", "_con")) %>%
  filter(n_good_ca >= n_cells, n_good_con >= n_cells) %>%
  pull(well_uid)
if (length(eligible) == 0) stop("No well has ", n_cells,
                                " good-paced cells in both modalities.")
set.seed(42)
target_well <- sample(sort(eligible), 1)
cat("  Selected well:", target_well, "\n")

# ---- Extract & normalise traces ----
extract_traces <- function(params, results, well, n, window_ms) {
  idx <- which(
    params$well_uid      == well &
    params$qc_pass       == TRUE &
    params$signal_class  == "good_paced"
  )
  cat(sprintf("  %d good-paced cells in %s\n", length(idx), well))
  if (length(idx) == 0) stop("No good-paced cells for well ", well)
  set.seed(42)
  picked <- if (length(idx) <= n) idx else sample(idx, n)

  bind_rows(lapply(seq_along(picked), function(j) {
    res <- results[[picked[j]]]
    if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
    tr <- filter(res$trace, time_ms >= window_ms[1], time_ms <= window_ms[2])
    if (nrow(tr) < 10) return(NULL)
    rng <- range(tr$signal_bc, na.rm = TRUE)
    if (diff(rng) < 1e-10) return(NULL)
    tr$signal_norm <- (tr$signal_bc - rng[1]) / diff(rng)
    tr$cell_label  <- paste0("Cell ", j)
    tr
  }))
}

cat("Calcium:\n");    ca_tr  <- extract_traces(ca_params,  ca_results,  target_well, n_cells, trace_window)
cat("Contraction:\n"); con_tr <- extract_traces(con_params, con_results, target_well, n_cells, trace_window)

# Fix facet order
ca_tr$cell_label  <- factor(ca_tr$cell_label,  levels = paste0("Cell ", seq_len(n_cells)))
con_tr$cell_label <- factor(con_tr$cell_label, levels = paste0("Cell ", seq_len(n_cells)))

# ---- Plot skeleton ----
make_grid_plot <- function(traces, line_color) {
  ggplot(traces, aes(x = time_ms / 1000, y = signal_norm)) +
    geom_line(color = line_color, linewidth = 0.55, alpha = 0.85) +
    facet_wrap(~ cell_label, ncol = n_col, scales = "fixed") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)), breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.08)), breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 7) +
    theme(
      strip.text       = element_text(size = 5.5, margin = margin(b = 1)),
      strip.background = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_line(linewidth = 0.3, color = "grey60"),
      panel.spacing    = unit(1.5, "mm"),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

p_ca  <- make_grid_plot(ca_tr,  "#1B6CA8")
p_con <- make_grid_plot(con_tr, "#C94040")

# ---- Save ----
trace_fig_dir <- file.path(pub_dir, "figure1_trace_illustration")
dir.create(trace_fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(trace_fig_dir, paste0("figure_trace_calcium_",     analysis_date, ".pdf")),
  p_ca,  width = 80, height = 40, units = "mm", device = "pdf"
)
ggsave(
  file.path(trace_fig_dir, paste0("figure_trace_contraction_", analysis_date, ".pdf")),
  p_con, width = 80, height = 40, units = "mm", device = "pdf"
)
cat("\nSaved to:", trace_fig_dir, "\n")
