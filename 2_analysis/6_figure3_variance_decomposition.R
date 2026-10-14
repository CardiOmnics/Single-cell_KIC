# =============================================================================
# 6_figure3_variance_decomposition.R
# Figure 3: Rationale for single-cell resolution
#
#   A  Variance decomposition — all measured parameters
#   B  Well averages | single-cell distributions for Ca2+ amplitude, Ca2+ CTD50,
#      contraction rise time and contraction decay time
#
# Sourced from run_pipeline.R after 5_plot_and_stats.R. Uses calcium_gp,
# contraction_gp, pub_dir and xlsx_path from that script.
#
# Output: 3_output/Figures/Publication/figure3_variance_decomposition/
#   figure3.pdf        combined (183 x 247 mm)
#   figure3_A.pdf      panel A
#   figure3_B.pdf      panel B
# Also appends a Variance_Decomp sheet to statistical_results_<date>.xlsx.
# =============================================================================

cat("\n========== 6_figure3_variance_decomposition.R ==========\n")

library(cowplot)

fig_out <- file.path(pub_dir, "figure3_variance_decomposition")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Panel A: Variance decomposition (ICC) — all parameters incl. contr_amp
# =============================================================================
cat("Building Panel A: variance decomposition...\n")

model_params_all <- list(
  list(col = "ca_amp",           label = "Ca amplitude",        data = "calcium"),
  list(col = "ca_CTD50_ms",      label = "Ca CTD50",            data = "calcium"),
  list(col = "ca_CTD90_ms",      label = "Ca CTD90",            data = "calcium"),
  list(col = "ca_rise_ms",       label = "Ca rise time",        data = "calcium"),
  list(col = "ca_decay_ms",      label = "Ca decay time",       data = "calcium"),
  list(col = "Cell.Area",        label = "Cell area",           data = "calcium"),
  list(col = "Nucleus.Area",     label = "Nucleus area (S1)",   data = "calcium"),
  list(col = "NC_ratio",         label = "N/C ratio",           data = "calcium"),
  list(col = "contr_amp",        label = "Contr. amplitude",    data = "contraction"),
  list(col = "contr_rise_ms",    label = "Contr. rise time",    data = "contraction"),
  list(col = "contr_decay_ms",   label = "Contr. decay time",   data = "contraction"),
  list(col = "Nucleus.Area",     label = "Nucleus area (S2)",   data = "contraction"),
  list(col = "Bodipy_intensity", label = "Bodipy intensity",    data = "contraction")
)

icc_results <- map_dfr(model_params_all, function(mp) {
  src    <- if (mp$data == "calcium") calcium_gp else contraction_gp
  df_icc <- src %>%
    select(well_uid, Condition, Cell_Line, all_of(mp$col)) %>%
    drop_na()
  if (nrow(df_icc) < 10) return(NULL)

  # Cell_Line has only 2 levels — a random effect would collapse to 0.
  # Treat it as fixed and derive its variance from the fixed-effect predictions.
  fit <- tryCatch(
    lmer(as.formula(paste0("`", mp$col,
         "` ~ Cell_Line + (1|Cell_Line:Condition) + (1|well_uid)")),
         data = df_icc),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  vc            <- as.data.frame(VarCorr(fit))
  var_condition <- vc$vcov[vc$grp == "Cell_Line:Condition"]
  var_well      <- vc$vcov[vc$grp == "well_uid"]
  var_residual  <- vc$vcov[vc$grp == "Residual"]
  # Variance between cell lines: take one predicted value per cell line so
  # both lines get equal weight, regardless of how many cells each has.
  # Use the population variance of the two fitted cell-line means here;
  # stats::var() would use the sample variance and double this component.
  cl_preds      <- tapply(predict(fit, re.form = ~0), df_icc$Cell_Line, mean)
  var_cell_line <- mean((cl_preds - mean(cl_preds))^2)
  var_total     <- var_cell_line + var_condition + var_well + var_residual

  tibble(
    parameter             = mp$label,
    pct_cell_line         = 100 * var_cell_line / var_total,
    pct_condition         = 100 * var_condition / var_total,
    pct_well              = 100 * var_well      / var_total,
    pct_residual          = 100 * var_residual  / var_total
  )
})

icc_long <- icc_results %>%
  pivot_longer(starts_with("pct_"), names_to = "source", values_to = "pct") %>%
  mutate(
    source = recode(source,
                    pct_cell_line = "Between cell lines",
                    pct_condition = "Between condition groups",
                    pct_well      = "Between wells",
                    pct_residual  = "Within-well (cell-to-cell)"),
    source    = factor(source, levels = c("Between cell lines",
                                          "Between condition groups",
                                          "Between wells",
                                          "Within-well (cell-to-cell)")),
    parameter = factor(parameter, levels = rev(icc_results$parameter))
  )

p_A <- ggplot(icc_long, aes(x = pct, y = parameter, fill = source)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Between cell lines"         = "#1A3A5C",
    "Between condition groups"   = "#2166AC",
    "Between wells"              = "#92C5DE",
    "Within-well (cell-to-cell)" = "#F4A582"
  )) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 102)) +
  labs(x = "% of total variance", y = NULL, fill = NULL) +
  journal_theme +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(3, "mm"),
    legend.key.width = unit(5, "mm")
  )

# =============================================================================
# Panels B–E: well-average vs single-cell
# facet_wrap(~Cell_Line) keeps both cell lines readable in separate facets.
# Each panel = ggarrange(p_well, p_cell, ncol = 2).
# =============================================================================
cat("Building Panels B-E: well-average vs single-cell...\n")

ca_well_avgs <- calcium_gp %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(across(c(ca_amp, ca_CTD50_ms), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

contr_well_avgs <- contraction_gp %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(across(c(contr_rise_ms, contr_decay_ms), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

figure_params <- list(
  list(col = "ca_amp",
       label = "Ca\u00b2\u207a amplitude (a.u.)",
       data  = "calcium"),
  list(col = "ca_CTD50_ms",
       label = "Ca\u00b2\u207a CTD50 (ms)",
       data  = "calcium"),
  list(col = "contr_rise_ms",
       label = "Contr. rise time (ms)",
       data  = "contraction"),
  list(col = "contr_decay_ms",
       label = "Contr. decay time (ms)",
       data  = "contraction")
)

# make_metric_plots returns list(p_well, p_cell) — two bare ggplot objects.
# All assembly happens in one flat plot_grid call, eliminating nested padding.
make_metric_plots <- function(par, show_col_titles = FALSE, show_x_labels = FALSE) {
  well_src <- if (par$data == "calcium") ca_well_avgs else contr_well_avgs
  cell_src <- if (par$data == "calcium") calcium_gp   else contraction_gp

  cell_vals   <- cell_src[[par$col]]
  cell_vals   <- cell_vals[is.finite(cell_vals)]
  y_lo        <- quantile(cell_vals, 0.005)
  y_hi        <- quantile(cell_vals, 0.995)
  y_range     <- y_hi - y_lo
  shared_ylim <- c(y_lo - y_range * 0.02, y_hi + y_range * 0.08)
  n_above     <- sum(cell_vals > y_hi, na.rm = TRUE)

  # Remove non-showing elements entirely (element_blank removes the space).
  # plot_grid rel_heights compensates so data panels stay equal in size.
  row_theme <- theme(
    strip.background = element_blank(),
    strip.text  = if (show_col_titles)
                    element_text(size = 9, face = "bold")
                  else element_blank(),
    plot.title  = if (show_col_titles)
                    element_text(size = 8, face = "plain", hjust = 0.5)
                  else element_blank(),
    axis.text.x  = if (show_x_labels)
                    element_text(angle = 35, hjust = 1, size = 7)
                   else element_blank(),
    axis.ticks.x = if (show_x_labels) element_line() else element_blank(),
    plot.margin      = margin(1, 2, 1, 2),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.position  = "none"
  )

  p_well <- ggplot(well_src,
                   aes(x = Condition, y = .data[[par$col]],
                       color = Condition)) +
    geom_point(size = 1.2, alpha = 0.8,
               position = position_jitter(width = 0.15, seed = 42)) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.5, linewidth = 0.4, color = "grey20") +
    coord_cartesian(ylim = shared_ylim) +
    facet_wrap(~Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(title = "Well averages", x = NULL, y = par$label) +
    journal_theme +
    row_theme

  p_cell <- ggplot(cell_src,
                   aes(x = Condition, y = .data[[par$col]],
                       fill = Condition, color = Condition)) +
    geom_rain(alpha     = 0.5,
              rain.side = "r",
              violin.args      = list(trim = TRUE, alpha = 0.45, color = NA),
              boxplot.args     = list(outlier.shape = NA, alpha = 0.6,
                                      color = "grey25"),
              boxplot.args.pos = list(width = 0.10),
              point.args       = list(alpha = 0.07, size = 0.2,
                                      color = "grey35")) +
    coord_cartesian(ylim = shared_ylim) +
    facet_wrap(~Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_fill_manual(values  = group_colors, drop = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(title = "Single-cell data", x = NULL, y = NULL) +
    journal_theme +
    row_theme

  if (n_above > 0) {
    p_cell <- p_cell +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("\u2191 %d above", n_above),
               hjust = 1.05, vjust = 1.4, size = 2.5, color = "grey50",
               fontface = "italic")
  }

  list(p_well = p_well, p_cell = p_cell)
}

n <- length(figure_params)
plot_pairs <- lapply(seq_along(figure_params), function(i) {
  make_metric_plots(figure_params[[i]],
                    show_col_titles = (i == 1),
                    show_x_labels   = (i == n))
})

# Flatten: well1, cell1, well2, cell2, ...
all_metric_plots <- unlist(lapply(plot_pairs, function(x) list(x$p_well, x$p_cell)),
                            recursive = FALSE)

# =============================================================================
# Assemble with cowplot::plot_grid — one flat call, no nesting, no padding
# =============================================================================
cat("Assembling figure...\n")

# rel_heights compensates for removed structural elements:
# row 1 = larger (strip + title added back), rows 2-3 = compact,
# row 4 = larger (angled x-axis labels added back).
# Tune the multipliers if data panels still look unequal.
panel_B <- plot_grid(plotlist    = all_metric_plots,
                     nrow        = 4,
                     ncol        = 2,
                     rel_widths  = c(1, 1),
                     rel_heights = c(1.25, 1, 1, 1.25),
                     align       = "v",
                     axis        = "l")

fig3 <- plot_grid(p_A, panel_B,
                  nrow        = 2,
                  rel_heights = c(1, 2.6),
                  labels      = c("A", "B"),
                  label_size  = 9,
                  label_fontface = "bold")

# =============================================================================
# Save
# =============================================================================
cat("Saving to:", fig_out, "\n")

ggsave(file.path(fig_out, "figure3.pdf"), fig3,
       width = 183, height = 247, units = "mm", device = "pdf")

ggsave(file.path(fig_out, "figure3_A.pdf"),
       p_A + theme(legend.position = "bottom"),
       width = 183, height = 80, units = "mm", device = "pdf")

ggsave(file.path(fig_out, "figure3_B.pdf"),
       panel_B,
       width = 183, height = 160, units = "mm", device = "pdf")

cat("  Saved: figure3.pdf, figure3_A.pdf, figure3_B.pdf\n")

# =============================================================================
# Append Variance_Decomp sheet to statistical_results.xlsx
# =============================================================================
# xlsx_path is the workbook written by 5_plot_and_stats.R in this run
if (exists("xlsx_path") && file.exists(xlsx_path)) {
  wb <- loadWorkbook(xlsx_path)

  var_decomp_tbl <- icc_results %>%
    mutate(across(starts_with("pct_"), ~ round(.x, 2))) %>%
    rename(
      `% Between cell lines`         = pct_cell_line,
      `% Between condition groups`   = pct_condition,
      `% Between wells`              = pct_well,
      `% Within-well (cell-to-cell)` = pct_residual
    )

  if ("Variance_Decomp" %in% names(wb)) removeWorksheet(wb, "Variance_Decomp")
  addWorksheet(wb, "Variance_Decomp")
  writeData(wb, "Variance_Decomp", relabel_cell_line_table(var_decomp_tbl))
  setColWidths(wb, "Variance_Decomp",
               cols = seq_along(var_decomp_tbl), widths = "auto")
  bold_style <- createStyle(textDecoration = "Bold")
  addStyle(wb, "Variance_Decomp", bold_style,
           rows = 1, cols = seq_along(var_decomp_tbl), gridExpand = TRUE)

  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("  Appended Variance_Decomp sheet to:", basename(xlsx_path), "\n")
} else {
  cat("  WARNING: statistical_results.xlsx from this run not found — Variance_Decomp not saved.\n")
}
