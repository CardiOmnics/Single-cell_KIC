# =============================================================================
# 5_plot_and_stats.R
# Merge all processed single-cell data, create publication-quality figures,
# correlation analyses, and mixed-effects models.
# Sourced from master pipeline.
# =============================================================================

if (!requireNamespace("ggcorrplot", quietly = TRUE)) install.packages("ggcorrplot")
if (!requireNamespace("ggrain",    quietly = TRUE)) install.packages("ggrain")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
if (!requireNamespace("cluster",   quietly = TRUE)) install.packages("cluster")
if (!requireNamespace("car",       quietly = TRUE)) install.packages("car")
if (!requireNamespace("moments",   quietly = TRUE)) install.packages("moments")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("cowplot",   quietly = TRUE)) install.packages("cowplot")
library(ggcorrplot)
library(ggrain)

cat("\n========== 5_plot_and_stats.R ==========\n")

# =============================================================================
# Section 1: Prepare session-aware datasets
# =============================================================================
# IMPORTANT: Calcium and Contraction+Bodipy were imaged in SEPARATE sessions.
# The plate shifted between sessions, so cell_uid does NOT match across channels.
#   Session 1: Calcium (Raw.Ci.on.Cm.API.F.t.)
#              + Cell.Area, Nucleus.Area, NC_ratio (segmentation correct here)
#   Session 2: Contraction (Calcium.Stretch.F.t.) + Bodipy (TMRM.API.F.t.)
#              + Nucleus.Area only (Cell.Area unreliable: drawn as circle around nucleus)
#
# Single-cell correlations are only valid WITHIN a session.
# Cross-session comparisons use WELL-LEVEL averages.
# =============================================================================
cat("\n--- Section 1: Preparing session-aware datasets ---\n")

# --- Session 1: Calcium + morphology (Cell.Area + Nucleus.Area both valid) ---
# Cell.Area, Nucleus.Area, and NC_ratio come directly from 2_process_calcium.R
# Main analysis uses good_paced cells only. Irregular cells are retained in
# ca_cell_params_all for the irregular-cell analysis (Section 5).
calcium_data <- ca_cell_params_all %>%
  filter(qc_pass, signal_class == "good_paced") %>%
  select(cell_uid, well_uid, Plate, Well, Condition, Cell_Line,
         ca_CTD50_ms = CTD50_ms, ca_CTD90_ms = CTD90_ms,
         ca_rise_ms = rise_time_ms, ca_decay_ms = decay_time_ms,
         ca_amp = mean_amp, ca_n_beats = n_beats, ca_ibi_cv = ibi_cv,
         ca_snr = snr_amp,
         Cell.Area, Nucleus.Area, NC_ratio)

# Back-compat alias — the rest of the script uses `calcium_gp` heavily.
calcium_gp <- calcium_data

# --- Session 2: Contraction + Nucleus.Area + Bodipy ---
# Nucleus.Area and Bodipy_intensity come directly from 1_process_contraction.R
# Cell.Area is NOT included: segmentation draws circle around nucleus — unreliable
# Good-paced cells only; irregular cells remain in contr_cell_params_all.
contraction_data <- contr_cell_params_all %>%
  filter(qc_pass, signal_class == "good_paced") %>%
  select(cell_uid, well_uid, Plate, Well, Condition, Cell_Line,
         contr_amp = mean_amp, contr_rise_ms = rise_time_ms,
         contr_decay_ms = decay_time_ms, contr_n_beats = n_beats,
         contr_ibi_cv = ibi_cv, contr_snr = snr_amp,
         Nucleus.Area, Bodipy_intensity)

# Back-compat alias — the rest of the script uses `contraction_gp` heavily.
contraction_gp <- contraction_data

# --- Well-level averages (for cross-session correlations) ---
ca_well_avg <- calcium_gp %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(
    n_ca_cells   = n(),
    ca_amp     = mean(ca_amp, na.rm = TRUE),
    ca_CTD50_ms  = mean(ca_CTD50_ms, na.rm = TRUE),
    ca_CTD90_ms  = mean(ca_CTD90_ms, na.rm = TRUE),
    ca_rise_ms   = mean(ca_rise_ms, na.rm = TRUE),
    ca_decay_ms  = mean(ca_decay_ms, na.rm = TRUE),
    Cell.Area    = mean(Cell.Area, na.rm = TRUE),
    Nucleus.Area.S1 = mean(Nucleus.Area, na.rm = TRUE),
    NC_ratio        = mean(NC_ratio, na.rm = TRUE),
    .groups = "drop"
  )

contr_well_avg <- contraction_gp %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(
    n_contr_cells    = n(),
    contr_amp        = mean(contr_amp, na.rm = TRUE),
    contr_rise_ms    = mean(contr_rise_ms, na.rm = TRUE),
    contr_decay_ms   = mean(contr_decay_ms, na.rm = TRUE),
    Nucleus.Area.S2  = mean(Nucleus.Area, na.rm = TRUE),
    Bodipy_intensity = mean(Bodipy_intensity, na.rm = TRUE),
    .groups = "drop"
  )

well_merged <- ca_well_avg %>%
  inner_join(contr_well_avg, by = c("well_uid", "Condition", "Cell_Line"))

# Ensure factor levels
ensure_factors <- function(df) {
  if (exists("group_levels")) df$Condition <- factor(df$Condition, levels = group_levels)
  df$Cell_Line <- factor(df$Cell_Line, levels = c("1.6", "1.6.24"))
  df
}
calcium_data        <- ensure_factors(calcium_data)
calcium_gp          <- ensure_factors(calcium_gp)
contraction_data        <- ensure_factors(contraction_data)
contraction_gp          <- ensure_factors(contraction_gp)
well_merged             <- ensure_factors(well_merged)

# Print summary
cat("\nSession-aware dataset summary:\n")
cat("  Session 1 - Calcium + morphology (Cell.Area + Nucleus.Area + NC_ratio):\n")
cat("    Good-paced cells:", nrow(calcium_data), "\n")
cat("  Session 2 - Contraction + Nucleus.Area + Bodipy:\n")
cat("    Good-paced cells:", nrow(contraction_data), "\n")
cat("  Cross-session well-level merge:", nrow(well_merged), "wells\n")
cat("\nGood-paced calcium cells per Condition x Cell Line:\n")
print(table(calcium_gp$Condition, calcium_gp$Cell_Line))
cat("\nGood-paced contraction cells per Condition x Cell Line:\n")
print(table(contraction_gp$Condition, contraction_gp$Cell_Line))

# Output directories
# Figures split into:
#   pub_dir  = main_dir/Figures/Publication  — multi-panel publication figures
#   expl_dir = main_dir/Figures/Exploratory  — supplementary figures, single
#              panels, diagnostics, trace galleries
if (!exists("pub_dir"))
  pub_dir  <- file.path(main_dir, "Figures", "Publication")
if (!exists("expl_dir"))
  expl_dir <- file.path(main_dir, "Figures", "Exploratory")
fig_dir   <- expl_dir                              # default for sub-figures
table_dir <- file.path(main_dir, "Tables")
dir.create(pub_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(expl_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

param_key <- function(spec) {
  paste(spec$data, spec$col, sep = "__")
}

# =============================================================================
# Section 2: Mixed-effects models (run BEFORE plots so p-values are available)
# =============================================================================
cat("\n--- Section 2: Mixed-effects models ---\n")

# Each parameter is assigned to the correct dataset based on imaging session.
# log = TRUE  -> primary model fitted on log10(outcome); residual skewness > 1
#                confirmed by Section 2b normality diagnostics.
# log = FALSE -> residual skewness < 1 on raw scale; raw model is adequate.
model_params <- list(
  # Calcium (Session 1 -> calcium_gp)
  list(col = "ca_amp",     label = "Calcium amplitude",      data = "calcium",      log = TRUE),
  list(col = "ca_CTD50_ms",  label = "Calcium CTD50",           data = "calcium",      log = TRUE),
  list(col = "ca_CTD90_ms",  label = "Calcium CTD90",           data = "calcium",      log = TRUE),
  list(col = "ca_rise_ms",   label = "Calcium rise time",       data = "calcium",      log = TRUE),
  list(col = "ca_decay_ms",  label = "Calcium decay time",      data = "calcium",      log = TRUE),
  # Morphology (Session 1 -> calcium_gp; Cell.Area + Nucleus.Area valid here)
  list(col = "Cell.Area",    label = "Cell area",               data = "calcium",      log = TRUE),
  list(col = "Nucleus.Area", label = "Nucleus area (S1)",       data = "calcium",      log = TRUE),
  list(col = "NC_ratio",     label = "N/C ratio",               data = "calcium",      log = FALSE),
  # Contraction (Session 2 -> contraction_gp)
  list(col = "contr_amp",      label = "Contraction amplitude", data = "contraction",  log = TRUE),
  list(col = "contr_rise_ms",  label = "Contraction rise time", data = "contraction",  log = TRUE),
  list(col = "contr_decay_ms", label = "Contraction decay time",data = "contraction",  log = TRUE),
  # Nucleus.Area Session 2 (valid; Cell.Area unreliable in contraction session)
  list(col = "Nucleus.Area", label = "Nucleus area (S2)",       data = "contraction",  log = TRUE),
  # Bodipy (Session 2 -> contraction_gp)
  list(col = "Bodipy_intensity", label = "Bodipy intensity",    data = "contraction",  log = TRUE)
)

emm_options(pbkrtest.limit = 50000)

all_anova   <- list()
all_posthoc <- list()
all_fixed   <- list()

for (mp in model_params) {
  use_log <- isTRUE(mp$log)
  mp_key <- param_key(mp)
  cat("  Fitting model for:", mp$label, if (use_log) "[log10]" else "[raw]", "\n")

  # Select correct dataset based on session
  src <- if (mp$data == "calcium") calcium_gp else contraction_gp

  df_model <- src %>%
    select(well_uid, Condition, Cell_Line, all_of(mp$col)) %>%
    drop_na() %>%
    mutate(
      Condition = factor(Condition, levels = group_levels),
      Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24"))
    )

  if (n_distinct(df_model$well_uid) < 3) {
    cat("    Skipped (fewer than 3 wells)\n")
    next
  }

  # Log-transform outcome if flagged (residual skewness > 1 on raw scale)
  if (use_log) {
    df_model <- df_model %>%
      filter(.data[[mp$col]] > 0) %>%
      mutate(.y = log10(.data[[mp$col]]))
    formula_str <- ".y ~ Condition * Cell_Line + (1 | well_uid)"
  } else {
    formula_str <- paste0("`", mp$col, "` ~ Condition * Cell_Line + (1 | well_uid)")
  }

  fit <- tryCatch(
    lmerTest::lmer(as.formula(formula_str), data = df_model),
    error = function(e) {
      cat("    Model failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) next

  # ANOVA
  aov <- as.data.frame(anova(fit))
  aov$parameter <- mp$label
  aov$term <- rownames(aov)
  aov$log_transformed <- use_log
  all_anova[[mp_key]] <- aov

  # Fixed effects
  fe <- broom.mixed::tidy(fit, effects = "fixed") %>%
    mutate(parameter = mp$label, log_transformed = use_log)
  all_fixed[[mp_key]] <- fe

  # Post-hoc: pairwise comparisons of Condition within each Cell_Line
  emm <- emmeans(fit, pairwise ~ Condition | Cell_Line, adjust = "tukey")
  ph <- as.data.frame(emm$contrasts) %>%
    mutate(parameter = mp$label, log_transformed = use_log)
  all_posthoc[[mp_key]] <- ph

  cat("    Done\n")
}

anova_results   <- bind_rows(all_anova)
fixed_results   <- bind_rows(all_fixed)
posthoc_results <- bind_rows(all_posthoc)

cat("  Models complete. Parameters fitted:", length(all_anova), "\n")

# =============================================================================
# Section 2b: Normality diagnostics — validate primary models
# =============================================================================
# Extracts residuals from the already-fitted primary models (Section 2) and
# checks whether the log10 transformation (applied to 10 of 13 parameters) has
# sufficiently normalised the residuals.
#
# Criterion: residual skewness.
# Shapiro-Wilk is NOT used — at N > 500 it rejects even trivially small
# deviations and its p-value is always ~0, making it uninformative.
# Threshold: |skewness| < 1 is considered acceptable for LME inference.
#
# Outputs: 3_output/QC/normality/
# =============================================================================
cat("\n--- Section 2b: Normality diagnostics (validating primary models) ---\n")

if (!requireNamespace("moments", quietly = TRUE)) install.packages("moments")
library(moments)

norm_dir <- file.path(main_dir, "QC", "normality")
dir.create(norm_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: one QQ plot panel
qq_panel <- function(resids, title, skew) {
  skew_ok <- abs(skew) < 1
  df <- data.frame(resid = resids[is.finite(resids)])
  ggplot(df, aes(sample = resid)) +
    stat_qq(size = 0.4, alpha = 0.3, colour = if (skew_ok) "#4393C3" else "#D6604D") +
    stat_qq_line(colour = "firebrick", linewidth = 0.6) +
    labs(title = sprintf("%s\nskew = %.2f%s", title, skew,
                         if (skew_ok) "" else "  *** HIGH ***"),
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_project() +
    theme(plot.title = element_text(size = 7.5,
                                    colour = if (skew_ok) "grey20" else "#B2182B"))
}

all_norm_diag <- list()

for (mp in model_params) {
  use_log  <- isTRUE(mp$log)
  mp_key   <- param_key(mp)
  src      <- if (mp$data == "calcium") calcium_gp else contraction_gp
  raw_vals <- src[[mp$col]]
  raw_vals <- raw_vals[is.finite(raw_vals)]
  if (length(raw_vals) < 10) next

  skew_outcome_raw <- tryCatch(moments::skewness(raw_vals), error = function(e) NA_real_)

  # Re-extract residuals from the primary fitted model
  # (re-fit here since we don't store the fit objects above)
  df_m <- src %>%
    select(well_uid, Condition, Cell_Line, all_of(mp$col)) %>%
    drop_na() %>%
    mutate(Condition = factor(Condition, levels = group_levels),
           Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24")))

  if (use_log) {
    df_m   <- df_m %>%
      filter(.data[[mp$col]] > 0) %>%
      mutate(.y = log10(.data[[mp$col]]))
    fmla   <- ".y ~ Condition * Cell_Line + (1 | well_uid)"
  } else {
    fmla   <- paste0("`", mp$col, "` ~ Condition * Cell_Line + (1 | well_uid)")
  }

  fit_check <- tryCatch(lmerTest::lmer(as.formula(fmla), data = df_m),
                        error = function(e) NULL)
  if (is.null(fit_check)) next

  resids     <- residuals(fit_check)
  skew_resid <- tryCatch(moments::skewness(resids[is.finite(resids)]), error = function(e) NA_real_)

  all_norm_diag[[mp_key]] <- tibble(
    parameter        = mp$label,
    model_scale      = if (use_log) "log10" else "raw",
    n_cells          = length(raw_vals),
    skew_outcome_raw = round(skew_outcome_raw, 2),
    skew_resid       = round(skew_resid, 2),
    resid_ok         = !is.na(skew_resid) && abs(skew_resid) < 1
  )

  # Save QQ plot
  p_qq <- qq_panel(resids,
                   sprintf("%s  [%s]", mp$label, if (use_log) "log10" else "raw"),
                   skew_resid)
  fname_qq <- file.path(norm_dir,
                        paste0("qq_", gsub("[^a-zA-Z0-9]", "_", mp$col), "_", analysis_date, ".pdf"))
  ggsave(fname_qq, p_qq, width = 5, height = 4, device = cairo_pdf)

  cat(sprintf("  %-30s [%s]  resid skew = %+.2f  %s\n",
              mp$label,
              if (use_log) "log10" else "raw  ",
              skew_resid,
              if (!is.na(skew_resid) && abs(skew_resid) < 1) "OK" else "*** HIGH ***"))
}

norm_summary <- bind_rows(all_norm_diag)

cat("\n  Summary (|skew| < 1 = acceptable for LME; SW p not shown — uninformative at N > 500):\n")
print(as.data.frame(norm_summary), row.names = FALSE)

n_ok      <- sum(norm_summary$resid_ok,  na.rm = TRUE)
n_flagged <- sum(!norm_summary$resid_ok, na.rm = TRUE)
cat(sprintf("\n  %d / %d parameters have acceptable residual skewness after transformation.\n",
            n_ok, nrow(norm_summary)))
if (n_flagged > 0) {
  cat("  Still flagged (consider reviewing):",
      paste(norm_summary$parameter[!norm_summary$resid_ok], collapse = ", "), "\n")
}

write.csv(relabel_cell_line_table(norm_summary),
          file.path(norm_dir, paste0("normality_summary_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved:", norm_dir, "\n  QQ plots: blue = |skew| < 1 (OK), red = |skew| >= 1\n\n")

# =============================================================================
# Section 3: Raincloud plot helper (with mixed-model p-values)
# =============================================================================
cat("\n--- Section 3: Raincloud helper + Figure 2A-B ---\n")

plot_raincloud <- function(data, y_col, y_label, title, log_scale = FALSE,
                           posthoc_df = NULL) {
  vals <- data[[y_col]]
  vals <- vals[is.finite(vals)]
  if (log_scale) vals <- vals[vals > 0]

  # Compute axis limits and padding in the scale-appropriate space
  lv_rng <- NA_real_  # only used for log scale step calculation below
  if (log_scale) {
    # Work in log10 space so padding/clipping are symmetric on the displayed axis
    lv      <- log10(vals)
    lv_lo   <- quantile(lv, 0.005, na.rm = TRUE)
    lv_hi   <- quantile(lv, 0.995, na.rm = TRUE)
    lv_rng  <- lv_hi - lv_lo
    y_lo          <- 10^lv_lo
    y_hi          <- 10^lv_hi
    y_n_pos       <- 10^(lv_lo - lv_rng * 0.12)
    y_bracket_top <- 10^(lv_hi + lv_rng * 0.38)
    ylim_use      <- c(10^(lv_lo - lv_rng * 0.20), y_bracket_top)
  } else {
    # Clip view to 0.5th–99.5th percentile
    y_lo    <- quantile(vals, 0.005, na.rm = TRUE)
    y_hi    <- quantile(vals, 0.995, na.rm = TRUE)
    y_range <- y_hi - y_lo
    # Space below for n labels, above for p-value brackets
    y_n_pos       <- y_lo - y_range * 0.12
    y_bracket_top <- y_hi + y_range * 0.38
    ylim_use      <- c(y_lo - y_range * 0.20, y_bracket_top)
  }
  n_above <- sum(vals > y_hi, na.rm = TRUE)
  n_below <- sum(vals < y_lo, na.rm = TRUE)

  # Per-group n labels
  n_labels <- data %>%
    filter(is.finite(.data[[y_col]])) %>%
    count(Condition, Cell_Line) %>%
    mutate(label = paste0("n=", n))

  p <- ggplot(data, aes(x = Condition, y = .data[[y_col]],
                        fill = Condition, color = Condition)) +
    geom_rain(
      alpha        = 0.5,
      rain.side    = "r",         # half-violin right, dots left
      violin.args  = list(trim = TRUE, alpha = 0.45, color = NA),
      boxplot.args = list(outlier.shape = NA, width = 0.10,
                          alpha = 0.6, color = "grey25"),
      point.args   = list(alpha = 0.07, size = 0.25, color = "grey35")
    ) +
    geom_text(data = n_labels,
              aes(x = Condition, y = y_n_pos, label = label),
              inherit.aes = FALSE,
              size = 2.8, color = "grey45", vjust = 0) +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_fill_manual(values  = group_colors, drop = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = title, x = NULL, y = y_label) +
    theme(legend.position = "none")

  if (log_scale) {
    p <- p +
      scale_y_log10(labels = scales::label_comma()) +
      coord_cartesian(ylim = c(max(ylim_use[1], .Machine$double.eps), ylim_use[2]))
  } else {
    p <- p + coord_cartesian(ylim = ylim_use)
  }

  # Annotate clipped outliers
  if (n_above > 0 || n_below > 0) {
    parts <- c()
    if (n_above > 0) parts <- c(parts, sprintf("\u2191 %d above", n_above))
    if (n_below > 0) parts <- c(parts, sprintf("\u2193 %d below", n_below))
    p <- p + annotate("text", x = Inf, y = Inf,
                      label = paste(parts, collapse = "  "),
                      hjust = 1.05, vjust = 1.4, size = 2.5,
                      color = "grey50", fontface = "italic")
  }

  # Mixed-model post-hoc p-value brackets
  if (!is.null(posthoc_df) && nrow(posthoc_df) > 0) {
    # step must be in data units; for log scale use a multiplicative factor
    step <- if (log_scale) y_hi * (10^(lv_rng * 0.09) - 1) else (y_hi - y_lo) * 0.09

    pv <- posthoc_df %>%
      mutate(
        group1  = sub(" - .*", "", contrast),
        group2  = sub(".* - ", "", contrast),
        p.adj   = p.value,
        p.label = case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01  ~ "**",
          p.adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      ) %>%
      filter(group1 == "MM" | group2 == "MM") %>%
      mutate(
        tmp_g1 = ifelse(group1 == "MM", group1, group2),
        tmp_g2 = ifelse(group1 == "MM", group2, group1),
        group1 = tmp_g1,
        group2 = tmp_g2
      ) %>%
      select(group1, group2, p.adj, p.label, Cell_Line)

    if (nrow(pv) > 0) {
      pv <- pv %>%
        group_by(Cell_Line) %>%
        mutate(y.position = y_hi + step * row_number()) %>%
        ungroup()

      p <- p + ggpubr::stat_pvalue_manual(
        pv, label = "p.label",
        step.increase = 0,
        tip.length = 0.01, size = 3.2
      )
    }
  }

  p
}

# =============================================================================
# Figure 2A-B: Troponin-I + BODIPY
#   A  Troponin-I barplot (n = 1 pooled measurement per condition/cell line)
#   B  BODIPY raincloud with mixed-model post-hoc brackets (vs MM)
# Panel C (immunostaining images) is added by hand.
# =============================================================================
cat("\n--- Figure: Troponin + Bodipy combined ---\n")

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)


# =============================================================================
# Panel A: Troponin-I barplot
# =============================================================================
troponin_raw <- read.csv(
  "../1_data/troponin.csv",
  header      = FALSE,
  skip        = 2,
  col.names   = c("short_name", "full_name", "troponin_ngL"),
  colClasses  = c("character", "character", "numeric")
)

troponin <- troponin_raw %>%
  mutate(
    Cell_Line = case_when(
      grepl("^1\\.6\\.24", short_name) ~ "1.6.24",
      TRUE                              ~ "1.6"
    ),
    cond_raw  = trimws(gsub("^1\\.6(\\.24)?\\s+", "", short_name)),
    Condition = case_when(
      cond_raw == "MM"    ~ "MM",
      cond_raw == "PA"    ~ "Palmitate",
      cond_raw == "OL"    ~ "Oleate",
      cond_raw == "PA+OL" ~ "Combination",
      TRUE                ~ cond_raw
    ),
    Condition = factor(Condition, levels = group_levels)
  )

p_troponin <- ggplot(troponin,
                     aes(x = Condition, y = troponin_ngL, fill = Condition)) +
  geom_col(width = 0.65, color = "grey25", linewidth = 0.3) +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
  labs(
    x       = NULL,
    y       = "Troponin-I (ng/L)",
    caption = "n=1; one pooled measurement from 6 wells per condition"
  ) +
  journal_theme +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1),
    plot.caption     = element_text(size = 6, color = "grey50", hjust = 0),
    legend.position  = "none"
  )

# =============================================================================
# Panel B: Bodipy raincloud (journal-themed)
# =============================================================================
# all_posthoc is keyed by param_key() = "<session>__<column>"
ph_bodipy_key <- "contraction__Bodipy_intensity"
if (!ph_bodipy_key %in% names(all_posthoc)) {
  stop("No post-hoc results for '", ph_bodipy_key, "' — Figure 2B needs them.")
}
ph_bodipy <- all_posthoc[[ph_bodipy_key]]

p_bodipy <- plot_raincloud(
  contraction_gp,
  y_col      = "Bodipy_intensity",
  y_label    = "Bodipy intensity (a.u.)",
  title      = "Bodipy lipid content",
  log_scale  = FALSE,
  posthoc_df = ph_bodipy
) +
  journal_theme +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    plot.title      = element_blank(),
    legend.position = "none"
  )

# =============================================================================
# Combine with patchwork
# =============================================================================
combined_fig <- (p_troponin / p_bodipy) +
  plot_layout(heights = c(1, 1.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

# =============================================================================
# Save
# =============================================================================
fig_out_dir <- file.path(pub_dir, "figure_troponin_bodipy")
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(fig_out_dir, paste0("figure_troponin_bodipy_", analysis_date, ".pdf")),
  combined_fig,
  width = 183, height = 130, units = "mm", device = "pdf"
)

ggsave(
  file.path(fig_out_dir, paste0("figure_troponin_bodipy_A_", analysis_date, ".pdf")),
  p_troponin,
  width = 183, height = 50, units = "mm", device = "pdf"
)

ggsave(
  file.path(fig_out_dir, paste0("figure_troponin_bodipy_B_", analysis_date, ".pdf")),
  p_bodipy,
  width = 183, height = 80, units = "mm", device = "pdf"
)

cat("  Saved: figure_troponin_bodipy/ (combined + A + B)\n")

# =============================================================================
# Section 4: Correlation analysis (session-aware)
# =============================================================================
cat("\n--- Section 4: Correlation analysis ---\n")

# --- Helper ---
cor_to_long <- function(df, cols, labels) {
  mat <- df %>% select(all_of(cols)) %>% drop_na()
  if (nrow(mat) < 10) return(NULL)
  colnames(mat) <- labels
  r <- cor(mat, method = "spearman")
  p <- cor_pmat(mat, method = "spearman")
  n <- nrow(mat)
  idx <- which(lower.tri(r), arr.ind = TRUE)
  tibble(
    var1 = factor(labels[idx[, 1]], levels = labels),
    var2 = factor(labels[idx[, 2]], levels = labels),
    rho = r[idx], pval = p[idx], sig = pval < 0.05, n = n
  )
}

plot_cor_heatmap <- function(cor_df, title, text_size = 2.5) {
  ggplot(cor_df, aes(x = var2, y = var1, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(sig, sprintf("%.2f", rho), "")), size = text_size) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         limits = c(-1, 1), name = "Spearman \u03C1") +
    theme_project() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    labs(title = title, x = NULL, y = NULL)
}

# --- 4a. Session 1 within-session correlations: Calcium + Morphology ---
cat("  Session 1 correlations (calcium + morphology, single-cell)...\n")

s1_cols <- c("ca_amp", "ca_CTD50_ms", "ca_CTD90_ms", "ca_rise_ms", "ca_decay_ms",
             "Cell.Area", "Nucleus.Area", "NC_ratio")
s1_labels <- c("Ca amplitude", "Ca CTD50", "Ca CTD90", "Ca rise", "Ca decay",
               "Cell area", "Nuc area", "N/C ratio")

cor_s1 <- cor_to_long(calcium_gp, s1_cols, s1_labels)
if (!is.null(cor_s1)) {
  p <- plot_cor_heatmap(cor_s1, paste0("Session 1: Calcium + Morphology (n=", cor_s1$n[1], " cells)"))
  ggsave(file.path(fig_dir, paste0("cor_session1_all_", analysis_date, ".pdf")),
         p, width = 8, height = 7, device = cairo_pdf)
  cat("  Saved: cor_session1_all\n")
}

# --- 4b. Session 2 within-session correlations: Contraction + Nucleus.Area + Bodipy ---
cat("  Session 2 correlations (contraction + Nucleus.Area + Bodipy, single-cell)...\n")

s2_cols <- c("contr_amp", "contr_rise_ms", "contr_decay_ms", "Nucleus.Area", "Bodipy_intensity")
s2_labels <- c("Contr amp", "Contr rise", "Contr decay", "Nuc area (S2)", "Bodipy")

cor_s2 <- cor_to_long(contraction_gp, s2_cols, s2_labels)
if (!is.null(cor_s2)) {
  p <- plot_cor_heatmap(cor_s2, paste0("Session 2: Contraction + Nucleus.Area + Bodipy (n=", cor_s2$n[1], " cells)"))
  ggsave(file.path(fig_dir, paste0("cor_session2_all_", analysis_date, ".pdf")),
         p, width = 7, height = 6, device = cairo_pdf)
  cat("  Saved: cor_session2_all\n")
}

# --- 4c. Well-level Ca x Contraction pairwise correlations (stats table) ---
# Point estimates + CIs + p-values for the specific cross-modality pairs
# reported in the paper. Uses `well_merged` (~34 wells with both modalities).
cat("  Well-level Ca x Contraction pairwise correlations (stats)...\n")

cross_pairs_to_test <- tibble::tribble(
  ~x_var,        ~y_var,
  "ca_decay_ms", "contr_amp",
  "ca_amp",      "contr_amp",
  "ca_CTD50_ms", "contr_amp",
  "ca_CTD90_ms", "contr_amp",
  "ca_rise_ms",  "contr_rise_ms",
  "ca_decay_ms", "contr_decay_ms",
  "ca_amp",      "contr_decay_ms"
)

well_cross_cor <- pmap_dfr(cross_pairs_to_test, function(x_var, y_var) {
  d <- well_merged %>% select(all_of(c(x_var, y_var))) %>% drop_na()
  if (nrow(d) < 5) return(NULL)
  ct <- cor.test(d[[x_var]], d[[y_var]], method = "spearman", exact = FALSE)
  tibble(
    x_var = x_var, y_var = y_var, n_wells = nrow(d),
    spearman_rho = round(unname(ct$estimate), 3),
    p_value      = signif(ct$p.value, 4)
  )
})
cat("\n  Well-level cross-modality correlations (Spearman):\n")
print(as.data.frame(well_cross_cor), row.names = FALSE)

# Faceted scatter plot with per-panel r/p labels
well_cross_long <- pmap_dfr(cross_pairs_to_test, function(x_var, y_var) {
  well_merged %>%
    transmute(well_uid, Condition, Cell_Line,
              x = .data[[x_var]], y = .data[[y_var]],
              facet = paste(y_var, "vs", x_var))
})
well_cross_labels <- well_cross_cor %>%
  mutate(facet = paste(y_var, "vs", x_var),
         label = sprintf("\u03C1 = %.2f, p = %.3g  (n=%d)",
                         spearman_rho, p_value, n_wells))

p_well_cross <- ggplot(well_cross_long,
                       aes(x = x, y = y,
                           color = Condition, shape = Cell_Line)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "grey30", linewidth = 0.5, alpha = 0.15) +
  facet_wrap(~ facet, scales = "free", ncol = 2) +
  geom_text(data = well_cross_labels,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.3, inherit.aes = FALSE, size = 2.8) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  theme_project() +
  labs(title = "Well-level calcium x contraction correlations",
       subtitle = paste0("Each point = one well (n = ", nrow(well_merged),
                         " wells with both modalities)"),
       x = NULL, y = NULL) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir,
                 paste0("fig_well_cross_correlations_",
                        analysis_date, ".pdf")),
       p_well_cross, width = 9, height = 10, device = cairo_pdf)
write.csv(relabel_cell_line_table(well_cross_cor),
          file.path(table_dir,
                    paste0("well_cross_correlations_",
                           analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: fig_well_cross_correlations + well_cross_correlations.csv\n")

# =============================================================================
# Section 5: Irregular cell analysis (session-aware)
# =============================================================================
cat("\n--- Section 5: Irregular cell analysis ---\n")

# Proportion of irregular cells per condition x cell line — computed per session.
# Sourced from the full *_cell_params_all tables (NOT calcium_data /
# contraction_data, which are good_paced-only). Denominator = all classifiable
# cells (qc_pass == TRUE, i.e. cells with a valid trace).
ca_source <- ca_cell_params_all %>%
  filter(qc_pass) %>%
  mutate(
    Condition = if (exists("group_levels"))
      factor(Condition, levels = group_levels) else Condition,
    Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24"))
  )
contr_source <- contr_cell_params_all %>%
  filter(qc_pass) %>%
  mutate(
    Condition = if (exists("group_levels"))
      factor(Condition, levels = group_levels) else Condition,
    Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24"))
  )

# Session 1: Calcium
ca_irregular_summary <- ca_source %>%
  group_by(Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = n_total - n_good_paced,
    pct_irregular = round(100 * n_irregular / n_total, 1),
    .groups = "drop"
  ) %>%
  mutate(session = "Calcium")

# Session 2: Contraction
contr_irregular_summary <- contr_source %>%
  group_by(Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = n_total - n_good_paced,
    pct_irregular = round(100 * n_irregular / n_total, 1),
    .groups = "drop"
  ) %>%
  mutate(session = "Contraction")

irregular_summary <- bind_rows(ca_irregular_summary, contr_irregular_summary)

cat("\n  Proportion of irregular cells per group:\n")
print(as.data.frame(irregular_summary), row.names = FALSE)

# --- Per-well dot plot: % irregular cells per well ---
# Each dot = one well. Wells excluded by well-level QC (< min_good_paced
# good_paced cells) are shown as red triangles; included wells as filled
# circles colored by condition.
cat("\n  Per-well irregular % dot plot...\n")

build_well_irreg <- function(df, qc_gate, session_label) {
  df %>%
    filter(qc_pass | well_uid %in% qc_gate$well_uid) %>%
    group_by(well_uid, Condition, Cell_Line) %>%
    summarise(
      n_total     = n(),
      n_irregular = sum(signal_class == "irregular"),
      pct_irregular = 100 * n_irregular / n_total,
      .groups = "drop"
    ) %>%
    left_join(qc_gate %>% select(well_uid, well_qc_pass),
              by = "well_uid") %>%
    mutate(
      well_qc_pass = replace_na(well_qc_pass, TRUE),
      session = session_label
    )
}

well_irreg_ca    <- build_well_irreg(
  ca_cell_params_all %>% mutate(
    Condition = if (exists("group_levels")) factor(Condition, levels = group_levels) else Condition,
    Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24"))),
  ca_well_qc_gate, "Calcium")
well_irreg_contr <- build_well_irreg(
  contr_cell_params_all %>% mutate(
    Condition = if (exists("group_levels")) factor(Condition, levels = group_levels) else Condition,
    Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24"))),
  contr_well_qc_gate, "Contraction")
well_irreg_all   <- bind_rows(well_irreg_ca, well_irreg_contr) %>%
  mutate(
    Condition = if (exists("group_levels"))
      factor(Condition, levels = group_levels) else factor(Condition),
    Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24"))
  )

p_well_irreg <- ggplot(well_irreg_all,
                        aes(x = Condition, y = pct_irregular)) +
  # Included wells: filled circles, colored by condition
  geom_jitter(data = well_irreg_all %>% filter(well_qc_pass),
              aes(fill = Condition),
              shape = 21, size = 2.5, width = 0.2, alpha = 0.8,
              colour = "grey30", stroke = 0.3) +
  # Excluded wells: red triangles

  geom_jitter(data = well_irreg_all %>% filter(!well_qc_pass),
              shape = 24, size = 3, width = 0.2,
              fill = "red", colour = "red", alpha = 0.8) +
  facet_grid(session ~ Cell_Line, labeller = labeller(Cell_Line = cell_line_fig_labeller)) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  theme_project() +
  labs(title = "Irregular cells per well",
       subtitle = "Red triangles = wells excluded from main analysis (failed well QC)",
       x = NULL, y = "Irregular cells (%)") +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, paste0("dotplot_irregular_per_well_", analysis_date, ".pdf")),
       p_well_irreg, width = 10, height = 7, device = cairo_pdf)
cat("  Saved: dotplot_irregular_per_well\n")

# --- 5a. Irregular proportion: logistic mixed-effects model ---
# Per cell line x modality, test whether the probability of being classified
# as irregular depends on Condition: glmer() logistic model with well as
# random intercept (cells are clustered within wells), likelihood-ratio test
# against the intercept-only model.
cat("\n  Irregular proportion tests (GLMM likelihood-ratio test)...\n")

prep_irreg <- function(df, modality) {
  df %>%
    filter(qc_pass, signal_class %in% c("good_paced", "irregular")) %>%
    mutate(
      is_irregular = as.integer(signal_class == "irregular"),
      Condition    = if (exists("group_levels"))
                       factor(Condition, levels = group_levels)
                     else factor(Condition),
      Cell_Line    = factor(Cell_Line, levels = c("1.6", "1.6.24")),
      modality     = modality
    )
}
ca_irr  <- prep_irreg(ca_cell_params_all,    "Calcium")
con_irr <- prep_irreg(contr_cell_params_all, "Contraction")

fit_glmer_per_line <- function(df, modality) {
  map_dfr(levels(df$Cell_Line), function(cl) {
    sub <- df %>% filter(Cell_Line == cl)
    if (n_distinct(sub$Condition) < 2 || nrow(sub) < 20) return(NULL)

    fit <- tryCatch(
      glmer(is_irregular ~ Condition + (1 | well_uid),
            data = sub, family = binomial,
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 2e5))),
      error = function(e) NULL, warning = function(w) NULL
    )
    if (!is.null(fit)) {
      fit0 <- tryCatch(
        glmer(is_irregular ~ 1 + (1 | well_uid),
              data = sub, family = binomial,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5))),
        error = function(e) NULL, warning = function(w) NULL
      )
      lrt <- if (!is.null(fit0)) anova(fit0, fit) else NULL
      lrt_chi <- if (!is.null(lrt)) lrt$Chisq[2] else NA_real_
      lrt_df  <- if (!is.null(lrt)) lrt$Df[2]    else NA_real_
      lrt_p   <- if (!is.null(lrt)) lrt$`Pr(>Chisq)`[2] else NA_real_
    } else {
      lrt_chi <- lrt_df <- lrt_p <- NA_real_
    }

    tibble(
      modality    = modality,
      Cell_Line   = cl,
      n_cells     = nrow(sub),
      n_wells     = n_distinct(sub$well_uid),
      n_irregular = sum(sub$is_irregular),
      glmer_ok    = !is.null(fit),
      LRT_chisq   = round(lrt_chi, 3),
      LRT_df      = lrt_df,
      LRT_p       = signif(lrt_p, 4)
    )
  })
}

irregular_tests <- bind_rows(
  fit_glmer_per_line(ca_irr,  "Calcium"),
  fit_glmer_per_line(con_irr, "Contraction")
)

# Fixed-effect estimates (log-odds) from the same GLMMs
irregular_glmm_fe <- map_dfr(
  list(list(df = ca_irr,  mod = "Calcium"),
       list(df = con_irr, mod = "Contraction")),
  function(x) {
    map_dfr(levels(x$df$Cell_Line), function(cl) {
      sub <- x$df %>% filter(Cell_Line == cl)
      if (n_distinct(sub$Condition) < 2 || nrow(sub) < 20) return(NULL)
      fit <- tryCatch(
        glmer(is_irregular ~ Condition + (1 | well_uid),
              data = sub, family = binomial,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5))),
        error = function(e) NULL, warning = function(w) NULL
      )
      if (is.null(fit)) return(NULL)
      broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
        mutate(modality = x$mod, Cell_Line = cl) %>%
        select(modality, Cell_Line, term, estimate, std.error,
               conf.low, conf.high, statistic, p.value)
    })
  }
)

cat("\n  Irregular proportion: global LRT (GLMM):\n")
print(as.data.frame(irregular_tests), row.names = FALSE)

write.csv(relabel_cell_line_table(irregular_tests),
          file.path(table_dir, paste0("irregular_tests_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(irregular_glmm_fe),
          file.path(table_dir, paste0("irregular_glmer_fixed_effects_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: irregular_tests + irregular_glmer_fixed_effects\n")

# =============================================================================
# Section 6: Single-cell analyses (Figures 4-6)
# =============================================================================
cat("\n--- Section 6: Single-cell analyses (Figures 4-6) ---\n")

# ---- 6a. Well-average vs mixed-model statistical comparison ----------------
cat("  6a: Well-average vs single-cell stats comparison...\n")

# Match the mixed-model scale/test as closely as possible on well averages:
# use log10 whenever the primary model does, and extract a Type III Condition
# test from the fixed-effects-only linear model.
wellavg_condition_anova <- function(df) {
  old_contrasts <- getOption("contrasts")
  on.exit(options(contrasts = old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  fit <- lm(y ~ Condition * Cell_Line, data = df)
  car::Anova(fit, type = 3)
}

# Run well-average linear models and compare them to the mixed-model results
wellavg_vs_mixed <- map_dfr(model_params, function(mp) {
  mp_key <- param_key(mp)
  # Select correct dataset based on session
  src <- if (mp$data == "calcium") calcium_gp else contraction_gp

  n_cells_model <- src %>%
    filter(is.finite(.data[[mp$col]]))
  if (isTRUE(mp$log)) {
    n_cells_model <- n_cells_model %>%
      filter(.data[[mp$col]] > 0)
  }

  # Well-average approach
  wa <- src %>%
    group_by(well_uid, Condition, Cell_Line) %>%
    summarise(y = mean(.data[[mp$col]], na.rm = TRUE), .groups = "drop") %>%
    drop_na()

  if (isTRUE(mp$log)) {
    wa <- wa %>%
      filter(y > 0) %>%
      mutate(y = log10(y))
  }

  if (nrow(wa) < 6) return(NULL)

  fit_wa <- tryCatch(
    wellavg_condition_anova(wa),
    error = function(e) NULL
  )
  if (is.null(fit_wa)) return(NULL)

  # Get the Condition p-value from the matched well-average linear model
  p_wa <- fit_wa["Condition", "Pr(>F)"]

  # Get the Condition p-value from mixed model (already computed)
  aov_mm <- all_anova[[mp_key]]
  if (is.null(aov_mm)) return(NULL)
  p_mm <- aov_mm[aov_mm$term == "Condition", "Pr(>F)"]

  tibble(
    parameter  = mp$label,
    session    = mp$data,
    n_wells    = nrow(wa),
    n_cells    = nrow(n_cells_model),
    wellavg_scale = ifelse(mp$log, "log10", "raw"),
    condition_test = "Adjusted for cell line + interaction",
    p_wellavg  = p_wa,
    p_mixed    = p_mm,
    sig_wellavg = ifelse(p_wa < 0.05, "*", "ns"),
    sig_mixed   = ifelse(p_mm < 0.05, "*", "ns"),
    gained     = ifelse(p_wa >= 0.05 & p_mm < 0.05, "YES", "")
  )
})

if (nrow(wellavg_vs_mixed) > 0) {
  cat("\n  Well-average linear model vs Single-cell mixed model:\n")
  print(as.data.frame(wellavg_vs_mixed), row.names = FALSE)

  n_gained <- sum(wellavg_vs_mixed$gained == "YES", na.rm = TRUE)
  subtitle_text <- if (n_gained > 0) {
    sprintf(
      "Matched outcome scale and Type III Condition test in both approaches; mixed model is significant for %d additional parameter%s",
      n_gained, ifelse(n_gained == 1, "", "s")
    )
  } else {
    "Matched outcome scale and Type III Condition test in both approaches"
  }

  # Power comparison plot (extended with cluster proportion rows) is assembled
  # after clustering runs — see section 6b continued below.
}

# =============================================================================
# Section 6b: Subpopulation clustering analysis
# =============================================================================
# Cluster cells by functional profiles to ask: does lipid overload create
# distinct subpopulations? This is a key argument for single-cell resolution.
# Run per session (calcium+morphology and contraction+Bodipy separately).
cat("\n--- Section 6b: Subpopulation clustering analysis ---\n")

# run_subpop_analysis() and plot_cluster_traces() are defined in 0_setup/cluster_helpers.R
source("../0_setup/cluster_helpers.R")

# --- Run for Session 1: Calcium + Morphology ---
cat("  Session 1: Calcium + Morphology\n")
s1_clust_cols <- c("ca_amp", "ca_CTD50_ms", "ca_CTD90_ms",
                   "ca_rise_ms", "ca_decay_ms",
                   "Cell.Area", "Nucleus.Area", "NC_ratio")
s1_clust_labels <- c("Ca amplitude", "CTD50", "CTD90",
                     "Ca rise", "Ca decay",
                     "Cell area", "Nuc area", "N/C ratio")

# --- Session 2 variables (defined early so k diagnostics can use them) ---
s2_clust_cols <- c("contr_amp", "contr_rise_ms", "contr_decay_ms",
                   "Nucleus.Area", "Bodipy_intensity")
s2_clust_labels <- c("Contr amp", "Contr rise", "Contr decay",
                     "Nuc area (S2)", "Bodipy")

# k per session — adjust here if one session shows severe imbalance
k_calcium    <- 3L
k_contraction <- 3L

# --- Cluster-number diagnostics (silhouette + elbow, k = 2..6) ---
# Runs k-means for k = 2..6 on the same z-scored matrices used by
# run_subpop_analysis() and reports within-cluster SS (elbow) and mean
# silhouette width. Silhouette is computed on a 5 k-cell subsample so the
# distance matrix stays tractable.
cat("\n  Cluster-number diagnostics (k = 2..6)...\n")
if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")
library(cluster)

k_diagnostics <- function(data, cols, session_label,
                          k_range = 2:6, sil_sample = 5000) {
  clust_data <- data %>% select(all_of(cols)) %>% drop_na()
  if (nrow(clust_data) < 50) return(NULL)
  clust_mat <- scale(as.matrix(clust_data))     # z-score standardization

  set.seed(42)
  sil_idx <- if (nrow(clust_mat) > sil_sample)
    sample(seq_len(nrow(clust_mat)), sil_sample)
  else seq_len(nrow(clust_mat))
  d_sub <- dist(clust_mat[sil_idx, , drop = FALSE])

  map_dfr(k_range, function(k) {
    set.seed(42)
    km  <- kmeans(clust_mat, centers = k, nstart = 50, iter.max = 200)
    sil <- silhouette(km$cluster[sil_idx], d_sub)
    tibble(session = session_label, k = k,
           wss = km$tot.withinss,
           mean_sil = mean(sil[, "sil_width"]),
           n_cells = nrow(clust_mat),
           n_sil_used = length(sil_idx))
  })
}

k_tbl <- bind_rows(
  k_diagnostics(calcium_gp,     s1_clust_cols, "Calcium"),
  k_diagnostics(contraction_gp, s2_clust_cols, "Contraction")
)

cat("\n  k diagnostics (WSS = elbow; mean_sil = silhouette):\n")
print(as.data.frame(k_tbl %>% mutate(
  wss = round(wss, 0), mean_sil = round(mean_sil, 3))),
  row.names = FALSE)

k_long <- k_tbl %>%
  select(session, k, wss, mean_sil) %>%
  pivot_longer(c(wss, mean_sil),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         wss      = "Within-cluster SS (elbow)",
                         mean_sil = "Mean silhouette width"))

p_k <- ggplot(k_long, aes(x = k, y = value)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "#D95F02") +
  facet_grid(metric ~ session, scales = "free_y") +
  scale_x_continuous(breaks = 2:6) +
  theme_project() +
  labs(title = "Cluster-number diagnostics",
       subtitle = "Dashed line = k = 3 used in manuscript",
       x = "Number of clusters (k)", y = NULL)
ggsave(file.path(fig_dir,
                 paste0("fig_k_diagnostics_", analysis_date, ".pdf")),
       p_k, width = 8, height = 6, device = cairo_pdf)
write.csv(relabel_cell_line_table(k_tbl),
          file.path(table_dir, paste0("k_diagnostics_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: fig_k_diagnostics + k_diagnostics.csv\n")

subpop_calcium <- run_subpop_analysis(
  calcium_gp, s1_clust_cols, s1_clust_labels,
  "Calcium", fig_dir, analysis_date,
  k = k_calcium,
  cluster_palette = cluster_colors_calcium,
  # Fixed numbering: 3 = slowest transient (highest CTD50), 1 = largest cells,
  # 2 = the remaining small, fast cluster (matches the manuscript)
  label_rule = list(list(label = 3L, var = "ca_CTD50_ms", pick = "max"),
                    list(label = 1L, var = "Cell.Area",   pick = "max"))
)
if (!is.null(ca_results_all)) {
  p_ca_traces <- plot_cluster_traces(
    subpop_calcium, ca_cell_params_all, ca_results_all,
    session_label = "Calcium", y_label = "Calcium signal (baseline-corrected)",
    n_per_cluster = 8, fig_dir = fig_dir, analysis_date = analysis_date,
    cluster_palette = cluster_colors_calcium
  )
} else {
  p_ca_traces <- NULL
  cat("  Skipping calcium cluster traces (step 2 was loaded from RDS)\n")
}

# --- Calcium cluster figure (combined) ---
if (!is.null(subpop_calcium)) {
  if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
  library(cowplot)

  journal_theme_clust <- journal_theme +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())

  p_A <- subpop_calcium$p_comp    + journal_theme_clust +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 7))
  p_B <- subpop_calcium$p_profile + journal_theme_clust +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 7),
          legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.key.height = unit(2.5, "mm"),
          legend.key.width  = unit(8, "mm"),
          legend.title      = element_text(size = 6),
          legend.text       = element_text(size = 6))
  p_C <- subpop_calcium$p_pca_clust + journal_theme_clust + theme(legend.position = "bottom")

  fig_clust_dir <- file.path(pub_dir, "figure_calcium_clusters")
  dir.create(fig_clust_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(p_ca_traces)) {
    p_D <- p_ca_traces + journal_theme_clust +
      theme(strip.text = element_text(size = 9), legend.position = "none")

    top_row <- plot_grid(p_A, p_B, ncol = 2, rel_widths = c(1.1, 1),
                         labels = c("A", "B"), label_size = 9, label_fontface = "bold")
    fig_calcium_clusters <- plot_grid(top_row, p_C, p_D, nrow = 3,
                                      rel_heights = c(1.4, 0.85, 0.9),
                                      labels = c("", "C", "D"),
                                      label_size = 9, label_fontface = "bold")
    fig_h_mm <- 210
  } else {
    top_row <- plot_grid(p_A, p_B, ncol = 2, rel_widths = c(1.1, 1),
                         labels = c("A", "B"), label_size = 9, label_fontface = "bold")
    fig_calcium_clusters <- plot_grid(top_row, p_C, nrow = 2, rel_heights = c(1.5, 0.85),
                                      labels = c("", "C"),
                                      label_size = 9, label_fontface = "bold")
    fig_h_mm <- 150
    cat("  Panel D omitted: no trace data (run pipeline with run_calcium = TRUE)\n")
  }

  ggsave(file.path(fig_clust_dir, "figure_calcium_clusters.pdf"),
         fig_calcium_clusters,
         width = 183, height = fig_h_mm, units = "mm", device = "pdf")
  ggsave(file.path(fig_clust_dir, "figure_calcium_clusters_A.pdf"),
         p_A, width = 89, height = 80, units = "mm", device = "pdf")
  ggsave(file.path(fig_clust_dir, "figure_calcium_clusters_B.pdf"),
         p_B, width = 89, height = 70, units = "mm", device = "pdf")
  ggsave(file.path(fig_clust_dir, "figure_calcium_clusters_C.pdf"),
         p_C, width = 183, height = 80, units = "mm", device = "pdf")
  if (!is.null(p_ca_traces)) {
    ggsave(file.path(fig_clust_dir, "figure_calcium_clusters_D.pdf"),
           p_D, width = 183, height = 70, units = "mm", device = "pdf")
  }
  cat("  Saved: figure_calcium_clusters (", fig_h_mm, "mm tall)\n")
}

# --- Run for Session 2: Contraction + Nucleus.Area + Bodipy ---
cat("\n  Session 2: Contraction + Nucleus.Area + Bodipy\n")
# s2_clust_cols / s2_clust_labels defined above, alongside s1, so the
# k diagnostics block can reference them.

subpop_contraction <- run_subpop_analysis(
  contraction_gp, s2_clust_cols, s2_clust_labels,
  "Contraction + Nucleus.Area + Bodipy", fig_dir, analysis_date,
  k = k_contraction,
  cluster_palette = cluster_colors_contraction,
  # Fixed numbering: 3 = lipid-loaded (highest BODIPY), 1 = highest contraction
  # amplitude, 2 = the remaining low-amplitude cluster (matches the manuscript)
  label_rule = list(list(label = 3L, var = "Bodipy_intensity", pick = "max"),
                    list(label = 1L, var = "contr_amp",        pick = "max"))
)
if (!is.null(contr_results_all)) {
  p_contr_traces <- plot_cluster_traces(
    subpop_contraction, contr_cell_params_all, contr_results_all,
    session_label = "Contraction + Bodipy",
    y_label = "Contraction (baseline-corrected)",
    n_per_cluster = 8, fig_dir = fig_dir, analysis_date = analysis_date,
    cluster_palette = cluster_colors_contraction
  )
} else {
  p_contr_traces <- NULL
  cat("  Skipping contraction cluster traces (step 1 was loaded from RDS)\n")
}

# --- Contraction cluster figure (combined) ---
if (!is.null(subpop_contraction)) {
  if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
  library(cowplot)

  journal_theme_clust <- journal_theme +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())

  p_A <- subpop_contraction$p_comp    + journal_theme_clust +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 7))
  p_B <- subpop_contraction$p_profile + journal_theme_clust +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 7),
          legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.key.height = unit(2.5, "mm"),
          legend.key.width  = unit(8, "mm"),
          legend.title      = element_text(size = 6),
          legend.text       = element_text(size = 6))
  p_C <- subpop_contraction$p_pca_clust + journal_theme_clust + theme(legend.position = "bottom")

  fig_clust_dir <- file.path(pub_dir, "figure_contraction_clusters")
  dir.create(fig_clust_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(p_contr_traces)) {
    p_D <- p_contr_traces + journal_theme_clust +
      theme(strip.text = element_text(size = 9), legend.position = "none")

    top_row <- plot_grid(p_A, p_B, ncol = 2, rel_widths = c(1.1, 1),
                         labels = c("A", "B"), label_size = 9, label_fontface = "bold")
    fig_contraction_clusters <- plot_grid(top_row, p_C, p_D, nrow = 3,
                                          rel_heights = c(1.4, 0.85, 0.9),
                                          labels = c("", "C", "D"),
                                          label_size = 9, label_fontface = "bold")
    fig_h_mm <- 210
  } else {
    top_row <- plot_grid(p_A, p_B, ncol = 2, rel_widths = c(1.1, 1),
                         labels = c("A", "B"), label_size = 9, label_fontface = "bold")
    fig_contraction_clusters <- plot_grid(top_row, p_C, nrow = 2, rel_heights = c(1.5, 0.85),
                                          labels = c("", "C"),
                                          label_size = 9, label_fontface = "bold")
    fig_h_mm <- 150
    cat("  Panel D omitted: no trace data (run pipeline with run_contraction = TRUE)\n")
  }

  ggsave(file.path(fig_clust_dir, "figure_contraction_clusters.pdf"),
         fig_contraction_clusters,
         width = 183, height = fig_h_mm, units = "mm", device = "pdf")
  ggsave(file.path(fig_clust_dir, "figure_contraction_clusters_A.pdf"),
         p_A, width = 89, height = 80, units = "mm", device = "pdf")
  ggsave(file.path(fig_clust_dir, "figure_contraction_clusters_B.pdf"),
         p_B, width = 89, height = 70, units = "mm", device = "pdf")
  ggsave(file.path(fig_clust_dir, "figure_contraction_clusters_C.pdf"),
         p_C, width = 183, height = 80, units = "mm", device = "pdf")
  if (!is.null(p_contr_traces)) {
    ggsave(file.path(fig_clust_dir, "figure_contraction_clusters_D.pdf"),
           p_D, width = 183, height = 70, units = "mm", device = "pdf")
  }
  cat("  Saved: figure_contraction_clusters (", fig_h_mm, "mm tall)\n")
}

# ---- 6b continued: Two-panel power comparison figure (Figure 6) ----
# Panel A: well-average vs single-cell mixed model p-values per parameter.
# Panel B: within-cluster condition effect per parameter — linear mixed model
#          (well as random intercept) on single-cell values within each cluster.
#          Falls back to KW on per-well means when the mixed model cannot be
#          estimated (too few wells or singular fit), flagged in the console.
if (nrow(wellavg_vs_mixed) > 0) {

  # Cluster palettes come from 0_setup/colors.R — one hue family per session,
  # matching the cluster composition / PCA / trace plots in cluster_helpers.R.
  # Panel B facets calcium and contraction but shares a single colour scale, so
  # the key has to be session-aware: keying on "Cluster 1" alone would paint
  # calcium cluster 1 and contraction cluster 1 identically and imply a
  # correspondence between two independent clusterings that does not exist.

  # --- Within-cluster condition effect per parameter ---
  # Primary: lmer(y ~ Condition + Cell_Line + (1|well_uid)) on single-cell
  #   values within cluster k — mirrors the main mixed model exactly.
  # Fallback: KW on per-well means when lmer is singular or fails to converge
  #   (can happen when a cluster has very few cells in some wells).
  within_cluster_lmer <- function(session_data, clust_result, session_params) {
    if (is.null(clust_result)) return(NULL)
    data_cl <- session_data %>%
      inner_join(clust_result$clust_data %>% select(cell_uid, cluster),
                 by = "cell_uid")
    map_dfr(session_params, function(mp) {
      map_dfr(sort(unique(data_cl$cluster)), function(k) {
        sub_k <- data_cl %>%
          filter(cluster == k, is.finite(.data[[mp$col]]))
        if (isTRUE(mp$log)) sub_k <- sub_k %>%
          filter(.data[[mp$col]] > 0) %>%
          mutate(!!mp$col := log10(.data[[mp$col]]))

        # Hard minimum for KW fallback: at least 2 conditions and 3 wells total
        n_wells_per_cond <- sub_k %>%
          group_by(Condition) %>%
          summarise(n_wells = n_distinct(well_uid), .groups = "drop")
        n_wells_total <- sum(n_wells_per_cond$n_wells)
        if (length(unique(sub_k$Condition)) < 2 || n_wells_total < 3) return(NULL)

        # Stricter gate for lmer: >=2 wells per condition and >=10 cells
        lmer_ok <- min(n_wells_per_cond$n_wells) >= 2 && nrow(sub_k) >= 10

        sub_k$y <- sub_k[[mp$col]]

        # Try mixed model if data are sufficient; fall back to KW otherwise
        fit <- if (lmer_ok) {
          tryCatch(
            lmer(y ~ Condition + Cell_Line + (1 | well_uid),
                 data = sub_k, REML = FALSE),
            error = function(e) NULL
          )
        } else NULL
        used_fallback <- FALSE
        p_val <- if (!is.null(fit) && !isSingular(fit)) {
          aov_cl <- tryCatch(
            as.data.frame(anova(fit, type = "III")),
            error = function(e) NULL
          )
          if (!is.null(aov_cl) && "Pr(>F)" %in% colnames(aov_cl) &&
              "Condition" %in% rownames(aov_cl)) {
            aov_cl["Condition", "Pr(>F)"]
          } else { NA_real_ }
        } else {
          used_fallback <- TRUE
          well_means <- sub_k %>%
            group_by(well_uid, Condition) %>%
            summarise(y = mean(y), .groups = "drop")
          kt <- tryCatch(kruskal.test(y ~ Condition, data = well_means),
                         error = function(e) NULL)
          if (!is.null(kt)) kt$p.value else NA_real_
        }

        if (is.na(p_val)) return(NULL)

        # --- Direction of the condition effect within this cluster ----------
        # P values alone do not say whether conditions move a parameter up or
        # down within a cluster (needed to judge e.g. "opposing effects across
        # subpopulations"). Export per-condition means on the original scale
        # and the % difference vs MM. Mixed model: estimated marginal means
        # averaged over cell line, Dunnett-adjusted contrasts vs MM.
        # KW fallback: median of per-well means, no contrast P values.
        back <- if (isTRUE(mp$log)) function(x) 10^x else identity
        dir_tbl <- if (!used_fallback) {
          em <- tryCatch(emmeans::emmeans(fit, ~ Condition), error = function(e) NULL)
          if (is.null(em)) NULL else {
            em_df <- as.data.frame(em)
            ctr <- as.data.frame(emmeans::contrast(em, method = "trt.vs.ctrl",
                                                   ref = which(levels(em_df$Condition) == "MM"),
                                                   adjust = "dunnett"))
            ctr$Condition <- sub(" - MM$", "", ctr$contrast)
            mm_mean <- back(em_df$emmean[em_df$Condition == "MM"])
            em_df %>%
              transmute(Condition = as.character(Condition), mean = back(emmean)) %>%
              left_join(ctr %>% select(Condition, p_vs_MM = p.value), by = "Condition") %>%
              mutate(pct_vs_MM = 100 * (mean / mm_mean - 1))
          }
        } else {
          wm <- sub_k %>% group_by(well_uid, Condition) %>%
            summarise(y = mean(y), .groups = "drop") %>%
            group_by(Condition) %>% summarise(mean = back(median(y)), .groups = "drop")
          mm_mean <- wm$mean[wm$Condition == "MM"]
          wm %>% mutate(Condition = as.character(Condition), p_vs_MM = NA_real_,
                        pct_vs_MM = if (length(mm_mean)) 100 * (mean / mm_mean - 1) else NA_real_)
        }
        if (!is.null(dir_tbl)) {
          wc_direction_rows[[length(wc_direction_rows) + 1]] <<- dir_tbl %>%
            mutate(parameter = mp$label, session = mp$data, cluster = k,
                   n_cells = nrow(sub_k), n_wells = n_distinct(sub_k$well_uid),
                   p_condition = p_val, method = ifelse(used_fallback, "KW fallback", "lmer"),
                   .before = 1)
        }
        if (used_fallback)
          cat(sprintf("    [fallback KW] %s / Cluster %s\n", mp$label, k))

        tibble(
          parameter     = mp$label,
          session       = mp$data,
          neg_log_p     = -log10(p_val + 1e-20),
          method        = paste0("Cluster ", k),
          used_fallback = used_fallback
        )
      })
    })
  }

  ca_params  <- Filter(function(mp) mp$data == "calcium",     model_params)
  con_params <- Filter(function(mp) mp$data == "contraction", model_params)

  cat("\n  Within-cluster condition effect (mixed model; fallback = KW):\n")
  wc_direction_rows <- list()
  wc_rows <- bind_rows(
    within_cluster_lmer(calcium_gp,     subpop_calcium,     ca_params),
    within_cluster_lmer(contraction_gp, subpop_contraction, con_params)
  )

  print(as.data.frame(wc_rows %>%
    mutate(p_approx = 10^(-neg_log_p),
           sig = case_when(
             p_approx < 0.001 ~ "***", p_approx < 0.01 ~ "**",
             p_approx < 0.05  ~ "*",   TRUE             ~ "ns"
           )) %>%
    select(parameter, method, sig, used_fallback)), row.names = FALSE)

  # Shared parameter order across both panels (calcium first, then contraction,
  # so each session forms a contiguous block on the y axis)
  param_levels <- rev(wellavg_vs_mixed$parameter)

  # Pretty session labels for the horizontal facet strips
  session_levels <- c("contraction", "calcium")
  session_labels <- c(contraction = "Contraction", calcium = "Calcium")

  # Journal theme for both panels; horizontal strip on top of each sub-row
  pwr_theme <- theme_bw() +
    theme(
      axis.text        = element_text(size = 7),
      axis.title       = element_text(size = 9),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 7),
      strip.text.y     = element_text(size = 8, face = "bold", angle = 0),
      strip.background = element_rect(fill = "grey92", colour = NA),
      axis.line        = element_line(linewidth = 0.5),
      axis.ticks       = element_line(linewidth = 0.5),
      panel.border     = element_rect(linewidth = 0.5),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      legend.position  = "bottom",
      legend.key.size  = unit(3.5, "mm"),
      panel.spacing.y  = unit(2, "mm")
    )

  # --- Panel A: well-average vs single-cell mixed model ---
  comp_models <- wellavg_vs_mixed %>%
    select(parameter, session, p_wellavg, p_mixed) %>%
    pivot_longer(c(p_wellavg, p_mixed), names_to = "m", values_to = "pvalue") %>%
    filter(!is.na(pvalue)) %>%
    transmute(
      parameter,
      session = factor(session, levels = session_levels),
      neg_log_p = -log10(pvalue + 1e-20),
      method    = ifelse(m == "p_wellavg",
                         "Well-average LM",
                         "Single-cell mixed model")
    ) %>%
    mutate(param_order = factor(parameter, levels = param_levels))

  p_A <- ggplot(comp_models,
                aes(x = neg_log_p, y = param_order,
                    color = method, shape = method)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    facet_grid(session ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(session = session_labels)) +
    scale_color_manual(
      values = c("Well-average LM"          = "#999999",
                 "Single-cell mixed model"   = "#000000"),
      name = NULL
    ) +
    scale_shape_manual(
      values = c("Well-average LM"          = 17,
                 "Single-cell mixed model"   = 16),
      name = NULL
    ) +
    pwr_theme +
    labs(
      x = expression(-log[10](italic(p))),
      y = NULL
    )

  # --- Panel B: within-cluster condition effect per parameter ---
  # Colour key = session + cluster number, so each facet draws from its own hue
  # family. Legend entries stay readable ("Calcium cluster 2").
  wc_plot_data <- wc_rows %>%
    mutate(param_order  = factor(parameter, levels = param_levels),
           session      = factor(session, levels = session_levels),
           cluster_num  = as.integer(gsub("Cluster (\\d+)", "\\1", method)),
           cluster_key  = paste0(session_labels[as.character(session)],
                                 " cluster ", cluster_num))

  clust_key_tbl <- wc_plot_data %>%
    distinct(session, cluster_num, cluster_key) %>%
    mutate(col = ifelse(as.character(session) == "calcium",
                        cluster_colors_calcium[cluster_num],
                        cluster_colors_contraction[cluster_num])) %>%
    arrange(session, cluster_num)

  clust_color_vals <- setNames(clust_key_tbl$col, clust_key_tbl$cluster_key)
  clust_shape_vals <- setNames(rep(18L, nrow(clust_key_tbl)),
                               clust_key_tbl$cluster_key)

  wc_plot_data <- wc_plot_data %>%
    mutate(cluster_key = factor(cluster_key, levels = clust_key_tbl$cluster_key))

  p_B <- ggplot(wc_plot_data,
                aes(x = neg_log_p, y = param_order,
                    color = cluster_key, shape = cluster_key)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    facet_grid(session ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(session = session_labels)) +
    scale_color_manual(values = clust_color_vals, name = NULL) +
    scale_shape_manual(values = clust_shape_vals, name = NULL) +
    pwr_theme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(
      x = expression(-log[10](italic(p))),
      y = NULL
    )

  # --- Assemble with cowplot ---
  if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
  library(cowplot)

  fig_pwr <- plot_grid(p_A, p_B,
                       ncol        = 2,
                       rel_widths  = c(1.4, 1),
                       labels      = c("A", "B"),
                       label_size  = 9,
                       label_fontface = "bold",
                       align       = "h",
                       axis        = "tb")

  fig_pwr_dir <- file.path(pub_dir, "figure_power_comparison")
  dir.create(fig_pwr_dir, recursive = TRUE, showWarnings = FALSE)

  n_rows <- length(param_levels)
  fig_h_mm <- 35 + n_rows * 6.5   # extra padding for the two facet strips

  ggsave(file.path(fig_pwr_dir, "figure_power_comparison.pdf"),
         fig_pwr,
         width = 183, height = fig_h_mm, units = "mm", device = "pdf")
  ggsave(file.path(fig_pwr_dir, "figure_power_comparison_A.pdf"),
         p_A,
         width = 105, height = fig_h_mm, units = "mm", device = "pdf")
  ggsave(file.path(fig_pwr_dir, "figure_power_comparison_B.pdf"),
         p_B + theme(axis.text.y = element_text(size = 7),
                     axis.ticks.y = element_line(linewidth = 0.5)),
         width = 90, height = fig_h_mm, units = "mm", device = "pdf")

  cat(sprintf("  Saved: figure_power_comparison (A=methods, B=clusters; rows split by Contraction/Calcium; %.0f mm tall)\n",
              fig_h_mm))
  if (any(wc_rows$used_fallback, na.rm = TRUE)) {
    cat("  Note: some clusters used KW fallback (see [fallback KW] lines above)\n")
  }
}

# =============================================================================
# Section 6c: Lipid-loaded cells vs their neighbours in oleate wells
# =============================================================================
# BODIPY and contraction come from the same cells (session 2), so within one
# well lipid-loaded cells can be compared with the unloaded cells next to them.
# Both groups share medium and imaging conditions, which removes between-well
# variation from the comparison.
#   - Lipid-loaded = BODIPY intensity > 2x the median of MM (control) cells of
#     the same cell line.
#   - Per oleate well: difference in mean log10(parameter), loaded minus other.
#   - Test: one-sample t test of the per-well differences against 0, so the
#     WELL is the unit of replication (cells are not treated as replicates).
#     Both lines pooled; the line difference is tested with a two-sample t test
#     (equal variances) on the same per-well differences.
# Percentages are back-transformed: 100 * (10^diff - 1).
cat("\n--- Section 6c: Lipid-loaded vs neighbour cells (oleate wells) ---\n")

ll_fold      <- 2
ll_condition <- "Oleate"
ll_params    <- c(contr_amp      = "Contraction amplitude",
                  contr_rise_ms  = "Contraction rise time",
                  contr_decay_ms = "Contraction decay time")

ll_threshold <- contraction_gp %>%
  filter(Condition == "MM") %>%
  group_by(Cell_Line) %>%
  summarise(bodipy_threshold = ll_fold * median(Bodipy_intensity, na.rm = TRUE),
            .groups = "drop")

ll_cells <- contraction_gp %>%
  filter(Condition == ll_condition, is.finite(Bodipy_intensity)) %>%
  inner_join(ll_threshold, by = "Cell_Line") %>%
  mutate(lipid_loaded = Bodipy_intensity > bodipy_threshold)

# Per-well differences; a well needs both loaded and unloaded cells
lipid_loaded_wells <- ll_cells %>%
  group_by(well_uid, Cell_Line) %>%
  filter(any(lipid_loaded), any(!lipid_loaded)) %>%
  summarise(
    n_cells          = n(),
    n_loaded         = sum(lipid_loaded),
    pct_loaded       = round(100 * mean(lipid_loaded), 1),
    bodipy_threshold = first(bodipy_threshold),
    across(all_of(names(ll_params)),
           ~ mean(log10(.x[lipid_loaded & .x > 0])) -
             mean(log10(.x[!lipid_loaded & .x > 0])),
           .names = "dlog10_{.col}"),
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("dlog10_"), ~ 100 * (10^.x - 1),
                .names = "pct_{sub('dlog10_', '', .col)}"))

lipid_loaded_test <- map_dfr(names(ll_params), function(p) {
  d  <- lipid_loaded_wells[[paste0("dlog10_", p)]]
  tt <- t.test(d)
  by_line <- lipid_loaded_wells %>%
    group_by(Cell_Line) %>%
    summarise(m = mean(.data[[paste0("dlog10_", p)]]), .groups = "drop")
  line_p <- if (n_distinct(lipid_loaded_wells$Cell_Line) == 2)
    t.test(as.formula(paste0("dlog10_", p, " ~ Cell_Line")),
           data = lipid_loaded_wells, var.equal = TRUE)$p.value else NA_real_
  tibble(
    parameter          = ll_params[[p]],
    n_wells            = length(d),
    n_wells_lower      = sum(d < 0),
    pct_diff           = 100 * (10^tt$estimate - 1),
    pct_ci_low         = 100 * (10^tt$conf.int[1] - 1),
    pct_ci_high        = 100 * (10^tt$conf.int[2] - 1),
    t                  = unname(tt$statistic),
    df                 = unname(tt$parameter),
    p_value            = tt$p.value,
    pct_diff_R403L     = 100 * (10^by_line$m[by_line$Cell_Line == "1.6"]    - 1),
    pct_diff_ISO       = 100 * (10^by_line$m[by_line$Cell_Line == "1.6.24"] - 1),
    p_line_difference  = line_p
  )
}) %>%
  mutate(across(starts_with("pct_"), ~ round(.x, 1)))

cat("\n  Lipid-loaded cells per oleate well:\n")
print(as.data.frame(lipid_loaded_wells %>%
  select(well_uid, Cell_Line, n_cells, n_loaded, pct_loaded)), row.names = FALSE)
cat("\n  Loaded vs other cells (one-sample t test on per-well differences):\n")
print(as.data.frame(lipid_loaded_test), row.names = FALSE)

# =============================================================================
# Section 7: Export results (session-aware)
# =============================================================================
cat("\n--- Section 7: Exporting results ---\n")

# Good-paced cells per session
write.csv(relabel_cell_line_table(calcium_gp),
          file.path(table_dir, paste0("calcium_good_paced_", analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(contraction_gp),
          file.path(table_dir, paste0("contraction_good_paced_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: calcium_good_paced + contraction_good_paced CSVs\n")

# Summary statistics per session
ca_params <- sapply(model_params[sapply(model_params, "[[", "data") == "calcium"], "[[", "col")
contr_params <- sapply(model_params[sapply(model_params, "[[", "data") == "contraction"], "[[", "col")

summary_stats_ca <- calcium_gp %>%
  group_by(Condition, Cell_Line) %>%
  summarise(
    n = n(),
    across(all_of(ca_params),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

summary_stats_contr <- contraction_gp %>%
  group_by(Condition, Cell_Line) %>%
  summarise(
    n = n(),
    across(all_of(contr_params),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

write.csv(relabel_cell_line_table(summary_stats_ca),
          file.path(table_dir, paste0("summary_statistics_calcium_", analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(summary_stats_contr),
          file.path(table_dir, paste0("summary_statistics_contraction_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: summary_statistics CSVs\n")

# =============================================================================
# Excel workbook: comprehensive statistical results
# =============================================================================

# Helper: write a sheet with bold headers and optional p-value color highlighting
write_sheet <- function(wb, sheet, df, p_cols = NULL) {
  df <- relabel_cell_line_table(df)
  writeData(wb, sheet, df)
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  # Bold header row
  header_style <- createStyle(textDecoration = "bold",
                              fgFill = "#D9E1F2", border = "Bottom")
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:n_cols, gridExpand = TRUE)

  # Freeze top row
  freezePane(wb, sheet, firstRow = TRUE)

  # Auto column widths
  setColWidths(wb, sheet, cols = 1:n_cols, widths = "auto")

  # Conditional formatting for p-value columns
  if (!is.null(p_cols) && n_rows > 0) {
    sig_style  <- createStyle(fgFill = "#FFEB9C", fontColour = "#9C5700")   # yellow: p < 0.05
    vsig_style <- createStyle(fgFill = "#C6EFCE", fontColour = "#276221")   # green:  p < 0.01
    ns_style   <- createStyle(fgFill = "#FFFFFF")

    for (col_name in p_cols) {
      if (!col_name %in% names(df)) next
      col_idx <- which(names(df) == col_name)
      pvals   <- df[[col_name]]
      for (r in seq_len(n_rows)) {
        pv <- pvals[r]
        if (!is.na(pv) && is.numeric(pv)) {
          style <- if (pv < 0.01) vsig_style else if (pv < 0.05) sig_style else ns_style
          addStyle(wb, sheet, style, rows = r + 1, cols = col_idx)
        }
      }
    }
  }
}

# Add significance stars to posthoc results
posthoc_formatted <- posthoc_results %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.10  ~ ".",
      TRUE            ~ "ns"
    )
  ) %>%
  select(parameter, Cell_Line, contrast, estimate, SE, df, t.ratio, p.value, sig)

# Build clean Results Summary: one row per parameter x cell line
# Columns: ANOVA F + p for each term, then significant posthoc pairs
anova_wide <- anova_results %>%
  select(parameter, term, `F value`, `Pr(>F)`) %>%
  mutate(term = gsub(":", "x", term)) %>%
  pivot_wider(
    names_from  = term,
    values_from = c(`F value`, `Pr(>F)`),
    names_glue  = "{term}_{.value}"
  )

# Cell counts per parameter
cell_counts <- map_dfr(model_params, function(mp) {
  src <- if (mp$data == "calcium") calcium_gp else contraction_gp
  df_sub <- src %>% select(well_uid, all_of(mp$col)) %>% drop_na()
  tibble(
    parameter = mp$label,
    session   = mp$data,
    n_cells   = nrow(df_sub),
    n_wells   = n_distinct(df_sub$well_uid)
  )
})

# Significant posthoc pairs per parameter x cell line (collapsed to text)
sig_pairs <- posthoc_formatted %>%
  filter(p.value < 0.05) %>%
  group_by(parameter, Cell_Line) %>%
  summarise(
    sig_pairs = paste(
      sprintf("%s (%s)", contrast, sig),
      collapse = "; "
    ),
    .groups = "drop"
  )

results_summary <- cell_counts %>%
  left_join(anova_wide, by = "parameter") %>%
  left_join(
    sig_pairs %>%
      pivot_wider(names_from = Cell_Line,
                  values_from = sig_pairs,
                  names_prefix = "sig_pairs_"),
    by = "parameter"
  ) %>%
  rename_with(~ gsub(" value", "", .), starts_with("Condition_F")) %>%
  arrange(session, parameter)

# Tidy column order for readability
summary_col_order <- c(
  "parameter", "session", "n_cells", "n_wells",
  intersect(c("Condition_F value", "Condition_Pr(>F)",
              "Cell_Line_F value", "Cell_Line_Pr(>F)",
              "ConditionxCell_Line_F value", "ConditionxCell_Line_Pr(>F)"),
            names(results_summary)),
  grep("^sig_pairs", names(results_summary), value = TRUE)
)
results_summary <- results_summary %>%
  select(any_of(summary_col_order), everything())

# Assemble workbook
# Note: Supplementary Table 1 is now written by run_pipeline.R to its own
# file (supplementary_table_1_<date>.xlsx) so it is generated regardless of
# which downstream steps are toggled on/off.
wb <- createWorkbook()

# Sheet order: summary first for easy reading
sheets <- list(
  Results_Summary        = list(df = results_summary,
                                p_cols = grep("Pr", names(results_summary), value = TRUE)),
  ANOVA                  = list(df = anova_results,      p_cols = "Pr(>F)"),
  Posthoc                = list(df = posthoc_formatted,  p_cols = "p.value"),
  Fixed_Effects          = list(df = fixed_results,      p_cols = "p.value"),
  Summary_Calcium        = list(df = summary_stats_ca,   p_cols = NULL),
  Summary_Contraction    = list(df = summary_stats_contr,p_cols = NULL),
  Irregular_Summary      = list(df = irregular_summary,  p_cols = NULL)
)

if (exists("wellavg_vs_mixed") && nrow(wellavg_vs_mixed) > 0) {
  sheets[["Power_Comparison"]] <- list(df = wellavg_vs_mixed,
                                        p_cols = c("p_wellavg", "p_mixed"))
}

# Direction of within-cluster condition effects (Fig. 6B)
if (exists("wc_direction_rows") && length(wc_direction_rows) > 0) {
  within_cluster_direction <- bind_rows(wc_direction_rows)
  sheets[["WithinCluster_Direction"]] <- list(df = within_cluster_direction,
                                              p_cols = c("p_condition", "p_vs_MM"))
}

# Lipid-loaded vs neighbour cells in oleate wells (Section 6c)
if (exists("lipid_loaded_test")) {
  sheets[["LipidLoaded_Test"]]    <- list(df = lipid_loaded_test,
                                          p_cols = c("p_value", "p_line_difference"))
  sheets[["LipidLoaded_PerWell"]] <- list(df = lipid_loaded_wells, p_cols = NULL)
}

# Subpopulation clustering results
if (!is.null(subpop_calcium)) {
  sheets[["Clust_Ca_Profiles"]] <- list(df = subpop_calcium$cluster_profiles, p_cols = NULL)
  sheets[["Clust_Ca_Composition"]] <- list(df = subpop_calcium$composition, p_cols = NULL)
  sheets[["Clust_Ca_WellKW"]] <- list(df = subpop_calcium$well_kw, p_cols = "p_value")
}
if (!is.null(subpop_contraction)) {
  sheets[["Clust_Contr_Profiles"]] <- list(df = subpop_contraction$cluster_profiles, p_cols = NULL)
  sheets[["Clust_Contr_Composition"]] <- list(df = subpop_contraction$composition, p_cols = NULL)
  sheets[["Clust_Contr_WellKW"]] <- list(df = subpop_contraction$well_kw, p_cols = "p_value")
}

# Extra diagnostic sheets
if (exists("k_tbl"))                sheets[["K_Diagnostics"]]      <- list(df = k_tbl,                    p_cols = NULL)
if (exists("irregular_tests"))      sheets[["Irregular_Tests"]]    <- list(df = irregular_tests,         p_cols = "LRT_p")
if (exists("irregular_glmm_fe"))    sheets[["Irregular_GLMM_FE"]]  <- list(df = irregular_glmm_fe,       p_cols = "p.value")
if (exists("well_cross_cor"))       sheets[["Well_Cross_Cor"]]     <- list(df = well_cross_cor,          p_cols = "p_value")

for (nm in names(sheets)) {
  addWorksheet(wb, nm)
  write_sheet(wb, nm, sheets[[nm]]$df, sheets[[nm]]$p_cols)
}

xlsx_path <- file.path(table_dir, paste0("statistical_results_", analysis_date, ".xlsx"))
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat("  Saved: statistical_results_", analysis_date, ".xlsx\n", sep = "")
cat("    Sheets: Results_Summary | ANOVA | Posthoc | Fixed_Effects |\n")
cat("            Summary_Calcium | Summary_Contraction | Irregular_Summary")
if (exists("wellavg_vs_mixed") && nrow(wellavg_vs_mixed) > 0) cat(" | Power_Comparison")
cat("\n")

# =============================================================================
# Completion
# =============================================================================
cat("\n========== 4_plot_and_stats.R complete ==========\n")
cat("Outputs:\n")
cat("  Figures (publication):", pub_dir, "\n")
cat("    - figure_troponin_bodipy/      Figure 2A-B\n")
cat("    - figure_calcium_clusters/     Figure 4\n")
cat("    - figure_contraction_clusters/ Figure 5\n")
cat("    - figure_power_comparison/     Figure 6\n")
cat("  Figures (supplementary):", expl_dir, "\n")
cat("    - dotplot_irregular_per_well   Supplementary Fig. 1\n")
cat("    - cor_session1_all / cor_session2_all  Supplementary Fig. 3A / 3B\n")
cat("    - fig_well_cross_correlations  Supplementary Fig. 4\n")
cat("  Tables:", table_dir, "\n")
cat("    - statistical_results_<date>.xlsx (ANOVA = Supp. Table 2, Posthoc = Supp. Table 3,\n")
cat("      K_Diagnostics = Supp. Table 4, plus all other reported tests)\n")
cat("    - calcium_good_paced / contraction_good_paced (per-cell data)\n")
cat("================================================\n")
