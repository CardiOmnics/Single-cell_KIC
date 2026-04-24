# =============================================================================
# 4_plot_and_stats.R
# Merge all processed single-cell data, create publication-quality figures,
# correlation analyses, and mixed-effects models.
# Sourced from master pipeline.
# =============================================================================

if (!requireNamespace("ggcorrplot", quietly = TRUE)) install.packages("ggcorrplot")
if (!requireNamespace("ggrain",    quietly = TRUE)) install.packages("ggrain")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
if (!requireNamespace("cluster",   quietly = TRUE)) install.packages("cluster")
if (!requireNamespace("car",       quietly = TRUE)) install.packages("car")
library(ggcorrplot)
library(ggrain)

cat("\n========== 4_plot_and_stats.R ==========\n")

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
# ca_cell_params_all for the arrhythmia % plot and the PCA signal-class plot.
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

# All classifiable calcium cells (good_paced + irregular), with the same
# renamed columns as calcium_data. Used ONLY by the PCA signal-class sub-plot
# to show where irregular cells sit in parameter space.
calcium_all_classes <- ca_cell_params_all %>%
  filter(qc_pass) %>%
  select(cell_uid, well_uid, Plate, Well, Condition, Cell_Line,
         ca_CTD50_ms = CTD50_ms, ca_CTD90_ms = CTD90_ms,
         ca_rise_ms = rise_time_ms, ca_decay_ms = decay_time_ms,
         ca_amp = mean_amp, ca_n_beats = n_beats, ca_ibi_cv = ibi_cv,
         ca_snr = snr_amp, ca_signal_class = signal_class,
         Cell.Area, Nucleus.Area, NC_ratio)

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

# All classifiable contraction cells, with renamed columns. Used ONLY by the
# PCA signal-class sub-plot.
contraction_all_classes <- contr_cell_params_all %>%
  filter(qc_pass) %>%
  select(cell_uid, well_uid, Plate, Well, Condition, Cell_Line,
         contr_amp = mean_amp, contr_rise_ms = rise_time_ms,
         contr_decay_ms = decay_time_ms, contr_n_beats = n_beats,
         contr_ibi_cv = ibi_cv, contr_snr = snr_amp,
         contr_signal_class = signal_class,
         Nucleus.Area, Bodipy_intensity)

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
calcium_all_classes <- ensure_factors(calcium_all_classes)
contraction_data        <- ensure_factors(contraction_data)
contraction_gp          <- ensure_factors(contraction_gp)
contraction_all_classes <- ensure_factors(contraction_all_classes)
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
fig_dir   <- file.path(main_dir, "Figures")
table_dir <- file.path(main_dir, "Tables")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
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
# Section 3: Publication violin plots (with mixed-model p-values)
# =============================================================================
cat("\n--- Section 3: Generating violin plots ---\n")

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

slugify_filename <- function(x) {
  x %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", ., perl = TRUE) %>%
    gsub("^_+|_+$", "", ., perl = TRUE)
}

# Define all violin plots — each with correct session source
violin_specs <- list(
  # Calcium (Session 1)
  list(col = "ca_amp",    label = "Amplitude (a.u.)", title = "Calcium amplitude",              log = FALSE, data = "calcium"),
  list(col = "ca_CTD50_ms", label = "CTD50 (ms)",      title = "Calcium transient duration 50%", log = FALSE, data = "calcium"),
  list(col = "ca_CTD90_ms", label = "CTD90 (ms)",      title = "Calcium transient duration 90%", log = FALSE, data = "calcium"),
  list(col = "ca_rise_ms",  label = "Rise time (ms)",  title = "Calcium rise time",              log = TRUE,  data = "calcium"),
  list(col = "ca_decay_ms", label = "Decay time (ms)", title = "Calcium decay time",             log = TRUE,  data = "calcium"),
  # Morphology (Session 1 — Cell.Area + Nucleus.Area + NC_ratio all valid)
  list(col = "Cell.Area",    label = "Area (\u00B5m\u00B2)", title = "Cell area",               log = TRUE,  data = "calcium"),
  list(col = "Nucleus.Area", label = "Area (\u00B5m\u00B2)", title = "Nucleus area (S1)",       log = TRUE,  data = "calcium"),
  list(col = "NC_ratio",     label = "N/C ratio",       title = "Nucleus/Cell area ratio", log = FALSE, data = "calcium"),
  # Contraction (Session 2)
  list(col = "contr_amp",      label = "Amplitude (a.u.)", title = "Contraction amplitude",  log = FALSE, data = "contraction"),
  list(col = "contr_rise_ms",  label = "Rise time (ms)",   title = "Contraction rise time",  log = TRUE,  data = "contraction"),
  list(col = "contr_decay_ms", label = "Decay time (ms)",  title = "Contraction decay time", log = TRUE,  data = "contraction"),
  # Nucleus.Area Session 2 (valid; Cell.Area unreliable in contraction session)
  list(col = "Nucleus.Area", label = "Area (\u00B5m\u00B2)", title = "Nucleus area (S2)", log = TRUE, data = "contraction"),
  # Bodipy (Session 2)
  list(col = "Bodipy_intensity", label = "Bodipy intensity (a.u.)", title = "Bodipy lipid content", log = FALSE, data = "contraction"),
  # Supplementary: pacing quality
  list(col = "ca_ibi_cv",    label = "IBI CV", title = "Calcium beat regularity",    log = FALSE, data = "calcium"),
  list(col = "contr_ibi_cv", label = "IBI CV", title = "Contraction beat regularity", log = FALSE, data = "contraction")
)

for (vs in violin_specs) {
  cat("  Plotting:", vs$title, "\n")
  vs_key <- param_key(vs)
  src <- if (vs$data == "calcium") calcium_gp else contraction_gp
  ph_sub <- if (vs_key %in% names(all_posthoc)) all_posthoc[[vs_key]] else NULL
  p <- plot_raincloud(src, vs$col, vs$label, vs$title,
                      log_scale = vs$log, posthoc_df = ph_sub)
  fname <- paste0("raincloud_", slugify_filename(vs$title), "_", analysis_date, ".pdf")
  ggsave(file.path(fig_dir, fname), p, width = 10, height = 5.5, device = cairo_pdf)
}

cat("  Raincloud plots saved to:", fig_dir, "\n")

# =============================================================================
# Combined Troponin + Bodipy figure
#   Panel A: Troponin-I barplot (n=1 pooled measurement per condition/cell line)
#   Panel B: Bodipy raincloud (journal-themed)
# Uses contraction_gp, all_posthoc, group_colors, group_levels, plot_raincloud(),
#   fig_dir, analysis_date, journal_theme — all defined earlier in this script.
# =============================================================================
cat("\n--- Figure: Troponin + Bodipy combined ---\n")

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ---- Panel A: Troponin-I barplot --------------------------------------------
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

# ---- Panel B: Bodipy raincloud (journal-themed) -----------------------------
ph_bodipy <- if ("Bodipy_intensity" %in% names(all_posthoc)) {
  all_posthoc[["Bodipy_intensity"]]
} else {
  NULL
}

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

# ---- Combine with patchwork + save ------------------------------------------
combined_fig <- (p_troponin / p_bodipy) +
  plot_layout(heights = c(1, 1.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

fig_out_dir <- file.path(fig_dir, "figure_troponin_bodipy")
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
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
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

# Overall
cor_s1 <- cor_to_long(calcium_gp, s1_cols, s1_labels)
if (!is.null(cor_s1)) {
  p <- plot_cor_heatmap(cor_s1, paste0("Session 1: Calcium + Morphology (n=", cor_s1$n[1], " cells)"))
  ggsave(file.path(fig_dir, paste0("cor_session1_all_", analysis_date, ".pdf")),
         p, width = 8, height = 7, device = cairo_pdf)
  cat("  Saved: cor_session1_all\n")
}

# By Condition
cor_s1_cond <- map_dfr(group_levels, function(cond) {
  out <- cor_to_long(calcium_gp %>% filter(Condition == cond), s1_cols, s1_labels)
  if (!is.null(out)) out$facet <- paste0(cond, " (n=", out$n[1], ")")
  out
})
if (nrow(cor_s1_cond) > 0) {
  p <- plot_cor_heatmap(cor_s1_cond, "Session 1 correlations by Condition", text_size = 1.8) +
    facet_wrap(~ facet, nrow = 1) +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 9, face = "bold"))
  ggsave(file.path(fig_dir, paste0("cor_session1_by_condition_", analysis_date, ".pdf")),
         p, width = 20, height = 6, device = cairo_pdf)
  cat("  Saved: cor_session1_by_condition\n")
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

# By Condition
cor_s2_cond <- map_dfr(group_levels, function(cond) {
  out <- cor_to_long(contraction_gp %>% filter(Condition == cond), s2_cols, s2_labels)
  if (!is.null(out)) out$facet <- paste0(cond, " (n=", out$n[1], ")")
  out
})
if (nrow(cor_s2_cond) > 0) {
  p <- plot_cor_heatmap(cor_s2_cond, "Session 2 correlations by Condition", text_size = 2) +
    facet_wrap(~ facet, nrow = 1) +
    theme(strip.text = element_text(size = 9, face = "bold"))
  ggsave(file.path(fig_dir, paste0("cor_session2_by_condition_", analysis_date, ".pdf")),
         p, width = 16, height = 5, device = cairo_pdf)
  cat("  Saved: cor_session2_by_condition\n")
}

# --- 4c. Cross-session correlations at WELL LEVEL ---
cat("  Cross-session correlations (well-level averages)...\n")

cross_cols <- c("ca_amp", "ca_CTD50_ms", "ca_CTD90_ms", "ca_rise_ms", "ca_decay_ms",
                "Cell.Area", "Nucleus.Area.S1", "NC_ratio",
                "contr_amp", "contr_rise_ms", "contr_decay_ms",
                "Nucleus.Area.S2", "Bodipy_intensity")
cross_labels <- c("Ca amplitude", "Ca CTD50", "Ca CTD90", "Ca rise", "Ca decay",
                  "Cell area", "Nuc area (S1)", "N/C ratio",
                  "Contr amp", "Contr rise", "Contr decay",
                  "Nuc area (S2)", "Bodipy")

cor_cross <- cor_to_long(well_merged, cross_cols, cross_labels)
if (!is.null(cor_cross)) {
  p <- plot_cor_heatmap(cor_cross,
    paste0("Cross-session correlations (WELL-LEVEL averages, n=", cor_cross$n[1], " wells)")) +
    theme(plot.title = element_text(size = 11))
  ggsave(file.path(fig_dir, paste0("cor_cross_session_well_", analysis_date, ".pdf")),
         p, width = 9, height = 8, device = cairo_pdf)
  cat("  Saved: cor_cross_session_well\n")
}

# By Condition (well-level)
cor_cross_cond <- map_dfr(group_levels, function(cond) {
  out <- cor_to_long(well_merged %>% filter(Condition == cond), cross_cols, cross_labels)
  if (!is.null(out)) out$facet <- paste0(cond, " (n=", out$n[1], " wells)")
  out
})
if (nrow(cor_cross_cond) > 0) {
  p <- plot_cor_heatmap(cor_cross_cond,
    "Cross-session correlations by Condition (well-level)", text_size = 1.8) +
    facet_wrap(~ facet, nrow = 1) +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 9, face = "bold"))
  ggsave(file.path(fig_dir, paste0("cor_cross_session_by_condition_", analysis_date, ".pdf")),
         p, width = 24, height = 7, device = cairo_pdf)
  cat("  Saved: cor_cross_session_by_condition\n")
}

# --- 4d. Well-level Ca x Contraction pairwise correlations (stats table) ---
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
# Section 5: Key scatter plots (session-aware)
# =============================================================================
cat("\n--- Section 5: Scatter plots ---\n")

# Single-cell scatter (within-session)
plot_scatter_cell <- function(data, x_col, y_col, x_label, y_label, title) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = Condition)) +
    geom_point(alpha = 0.15, size = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = title, x = x_label, y = y_label)
}

# Well-level scatter (cross-session)
plot_scatter_well <- function(data, x_col, y_col, x_label, y_label, title) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = Condition)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = paste0(title, " (well averages)"),
         x = x_label, y = y_label)
}

# Within-session scatter plots (single-cell)
cat("  Within-session scatter plots (single-cell)...\n")
within_scatters <- list(
  # Session 2: Bodipy vs contraction
  list(data = "contraction", x = "Bodipy_intensity", y = "contr_amp",
       xl = "Bodipy intensity (a.u.)", yl = "Amplitude (a.u.)",
       title = "Bodipy vs contractility"),
  list(data = "contraction", x = "Bodipy_intensity", y = "contr_decay_ms",
       xl = "Bodipy intensity (a.u.)", yl = "Decay time (ms)",
       title = "Bodipy vs contraction decay"),
  # Session 1: Cell size vs calcium
  list(data = "calcium", x = "Cell.Area", y = "ca_amp",
       xl = "Cell area (\u00B5m\u00B2)", yl = "Ca amplitude (a.u.)",
       title = "Cell size vs calcium amplitude"),
  list(data = "calcium", x = "Cell.Area", y = "ca_CTD90_ms",
       xl = "Cell area (\u00B5m\u00B2)", yl = "CTD90 (ms)",
       title = "Cell size vs calcium reuptake")
)

for (ss in within_scatters) {
  cat("    ", ss$title, "\n")
  src <- if (ss$data == "calcium") calcium_gp else contraction_gp
  p <- plot_scatter_cell(src, ss$x, ss$y, ss$xl, ss$yl, ss$title)
  fname <- paste0("scatter_sc_", ss$x, "_vs_", ss$y, "_", analysis_date, ".pdf")
  ggsave(file.path(fig_dir, fname), p, width = 10, height = 5, device = cairo_pdf)
}

# Cross-session: Bodipy vs calcium, Contraction vs calcium (well-level)
cat("  Cross-session scatter plots (well averages)...\n")
cross_scatters <- list(
  list(x = "Bodipy_intensity", y = "ca_CTD90_ms",
       xl = "Bodipy intensity (a.u.)", yl = "CTD90 (ms)",
       title = "Bodipy vs calcium reuptake"),
  list(x = "Bodipy_intensity", y = "ca_amp",
       xl = "Bodipy intensity (a.u.)", yl = "Ca amplitude (a.u.)",
       title = "Bodipy vs calcium amplitude"),
  list(x = "contr_amp", y = "ca_amp",
       xl = "Contraction amplitude (a.u.)", yl = "Ca amplitude (a.u.)",
       title = "Contraction vs calcium amplitude"),
  list(x = "contr_amp", y = "ca_CTD90_ms",
       xl = "Contraction amplitude (a.u.)", yl = "CTD90 (ms)",
       title = "Contraction vs calcium reuptake")
)

for (ss in cross_scatters) {
  cat("    ", ss$title, "(well-level)\n")
  p <- plot_scatter_well(well_merged, ss$x, ss$y, ss$xl, ss$yl, ss$title)
  fname <- paste0("scatter_well_", ss$x, "_vs_", ss$y, "_", analysis_date, ".pdf")
  ggsave(file.path(fig_dir, fname), p, width = 10, height = 5, device = cairo_pdf)
}

cat("  Scatter plots saved to:", fig_dir, "\n")

# =============================================================================
# Section 6: Irregular cell analysis (session-aware)
# =============================================================================
cat("\n--- Section 6: Irregular cell analysis ---\n")

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

# Bar plot: % irregular cells per condition, faceted by session + cell line
p_irreg <- ggplot(irregular_summary,
                  aes(x = Condition, y = pct_irregular, fill = Condition)) +
  geom_col(alpha = 0.8, color = "grey30", width = 0.7) +
  geom_text(aes(label = paste0(pct_irregular, "%\n(", n_irregular, "/", n_total, ")")),
            vjust = -0.3, size = 2.8) +
  facet_grid(session ~ Cell_Line, labeller = labeller(Cell_Line = cell_line_fig_labeller)) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  theme_project() +
  labs(title = "Proportion of irregular/arrhythmic cells",
       subtitle = "Cells passing basic QC but without clean paced rhythm",
       x = NULL, y = "Irregular cells (%)") +
  theme(legend.position = "none") +
  ylim(0, max(irregular_summary$pct_irregular, na.rm = TRUE) * 1.25)
ggsave(file.path(fig_dir, paste0("barplot_irregular_pct_", analysis_date, ".pdf")),
       p_irreg, width = 10, height = 7, device = cairo_pdf)
cat("  Saved: barplot_irregular_pct\n")

# Save irregular cells data (per session)
ca_irregular_cells    <- ca_source    %>% filter(signal_class == "irregular")
contr_irregular_cells <- contr_source %>% filter(signal_class == "irregular")
write.csv(relabel_cell_line_table(ca_irregular_cells),
          file.path(table_dir, paste0("irregular_cells_calcium_", analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(contr_irregular_cells),
          file.path(table_dir, paste0("irregular_cells_contraction_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: irregular_cells_calcium + irregular_cells_contraction\n")

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

# --- 6a. Per-well cell-count summary (post-QC) ---
# Min / max / median / IQR of cells per well for each modality.
cat("\n  Per-well cell-count summary (post-QC)...\n")

well_counts_summary <- function(df, modality) {
  per_well <- df %>%
    filter(qc_pass) %>%
    count(well_uid, Condition, Cell_Line, name = "n_cells")
  per_well_gp <- df %>%
    filter(qc_pass, signal_class == "good_paced") %>%
    count(well_uid, Condition, Cell_Line, name = "n_good_paced") %>%
    right_join(per_well, by = c("well_uid", "Condition", "Cell_Line")) %>%
    mutate(n_good_paced = replace_na(n_good_paced, 0L))

  summ <- per_well_gp %>%
    summarise(
      modality     = modality,
      n_wells      = n(),
      min_cells    = min(n_cells),
      q25_cells    = quantile(n_cells, 0.25),
      median_cells = median(n_cells),
      q75_cells    = quantile(n_cells, 0.75),
      max_cells    = max(n_cells),
      iqr_cells    = IQR(n_cells),
      mean_cells   = round(mean(n_cells), 1),
      sd_cells     = round(sd(n_cells), 1),
      min_gp       = min(n_good_paced),
      median_gp    = median(n_good_paced),
      max_gp       = max(n_good_paced),
      iqr_gp       = IQR(n_good_paced)
    )

  list(summary = summ,
       per_well = per_well_gp %>% mutate(modality = modality))
}

wc_ca  <- well_counts_summary(ca_cell_params_all,    "Calcium")
wc_con <- well_counts_summary(contr_cell_params_all, "Contraction")

well_counts_overall <- bind_rows(wc_ca$summary, wc_con$summary)
well_counts_detail  <- bind_rows(wc_ca$per_well, wc_con$per_well)

cat("\n  Per-well cell counts after QC:\n")
print(as.data.frame(well_counts_overall), row.names = FALSE)

write.csv(relabel_cell_line_table(well_counts_overall),
          file.path(table_dir, paste0("per_well_cellcounts_summary_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(well_counts_detail),
          file.path(table_dir, paste0("per_well_cellcounts_detail_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: per_well_cellcounts_summary + per_well_cellcounts_detail\n")

# --- 6b. Irregular proportion: logistic mixed-effects + chi-square ---
# Per cell line x modality, test whether the probability of being classified
# as irregular depends on Condition. Primary: glmer() logistic model with
# well as random intercept (cells are clustered within wells). Fallback:
# chi-square on the raw contingency table (ignores well clustering).
cat("\n  Irregular proportion tests (GLMM + chi-square)...\n")

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

    ct  <- table(sub$Condition, sub$is_irregular)
    chi <- tryCatch(chisq.test(ct), error = function(e) NULL)

    tibble(
      modality    = modality,
      Cell_Line   = cl,
      n_cells     = nrow(sub),
      n_wells     = n_distinct(sub$well_uid),
      n_irregular = sum(sub$is_irregular),
      glmer_ok    = !is.null(fit),
      LRT_chisq   = round(lrt_chi, 3),
      LRT_df      = lrt_df,
      LRT_p       = signif(lrt_p, 4),
      chisq       = if (!is.null(chi)) round(chi$statistic, 3) else NA_real_,
      chisq_df    = if (!is.null(chi)) chi$parameter else NA_real_,
      chisq_p     = if (!is.null(chi)) signif(chi$p.value, 4) else NA_real_
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

cat("\n  Irregular proportion: global LRT (GLMM) + chi-square fallback:\n")
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

# --- 6c. Per-well breakdown of 1.6.24 control irregular calcium rate ---
# Is the elevated irregular calcium rate in 1.6.24 controls driven by
# specific wells, or spread evenly across replicates?
cat("\n  1.6.24 control: irregular calcium per well...\n")

ctrl_levels <- grep("control|ctrl|Control|CTRL",
                    levels(ca_irr$Condition), value = TRUE)
if (length(ctrl_levels) == 0) ctrl_levels <- levels(ca_irr$Condition)[1]
cat("  Control level(s):", paste(ctrl_levels, collapse = ", "), "\n")

ca_1624_ctrl_by_well <- ca_irr %>%
  filter(Cell_Line == "1.6.24", Condition %in% ctrl_levels) %>%
  group_by(well_uid, Plate, Well, Condition) %>%
  summarise(
    n_cells       = n(),
    n_irregular   = sum(is_irregular),
    pct_irregular = round(100 * n_irregular / n_cells, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_irregular))

cat("\n  Per-well irregular % in 1.6.24 controls:\n")
print(as.data.frame(ca_1624_ctrl_by_well), row.names = FALSE)
if (nrow(ca_1624_ctrl_by_well) > 0) {
  cat(sprintf("\n  Median = %.1f%%, IQR = %.1f%%, range = %.1f-%.1f%%\n",
              median(ca_1624_ctrl_by_well$pct_irregular),
              IQR(ca_1624_ctrl_by_well$pct_irregular),
              min(ca_1624_ctrl_by_well$pct_irregular),
              max(ca_1624_ctrl_by_well$pct_irregular)))

  p_1624 <- ggplot(ca_1624_ctrl_by_well,
                   aes(x = reorder(well_uid, pct_irregular),
                       y = pct_irregular)) +
    geom_col(fill = "#2166AC", alpha = 0.85) +
    geom_text(aes(label = paste0(n_irregular, "/", n_cells)),
              hjust = -0.1, size = 2.8) +
    coord_flip() +
    ylim(0, max(ca_1624_ctrl_by_well$pct_irregular) * 1.25) +
    theme_project() +
    labs(title = "1.6.24 control: irregular calcium cells per well",
         subtitle = "Is the elevated rate driven by specific wells?",
         x = NULL, y = "Irregular cells (%)")
  ggsave(file.path(fig_dir,
                   paste0("fig_1624_ctrl_irregular_by_well_",
                          analysis_date, ".pdf")),
         p_1624,
         width = 7,
         height = 2 + 0.3 * nrow(ca_1624_ctrl_by_well),
         device = cairo_pdf)
  write.csv(relabel_cell_line_table(ca_1624_ctrl_by_well),
            file.path(table_dir,
                      paste0("irregular_1624_ctrl_by_well_",
                             analysis_date, ".csv")),
            row.names = FALSE)
  cat("  Saved: fig_1624_ctrl_irregular_by_well + csv\n")
}

# =============================================================================
# Section 7: Single-cell manuscript figures
# =============================================================================
cat("\n--- Section 7: Single-cell manuscript figures ---\n")

# ---- 7a. "Well-average illusion" figure ------------------------------------
# Show well averages (like a traditional paper) next to single-cell distributions
cat("  7a: Well-average illusion figure...\n")

# Compute well averages per session
ca_well_avgs <- calcium_gp %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(
    n_cells = n(),
    across(c(ca_amp, ca_CTD50_ms, ca_CTD90_ms, ca_rise_ms, ca_decay_ms,
             Cell.Area, Nucleus.Area, NC_ratio),
           ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

contr_well_avgs <- contraction_gp %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(
    n_cells = n(),
    across(c(contr_amp, contr_rise_ms, contr_decay_ms, Nucleus.Area, Bodipy_intensity),
           ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Pick 4 representative parameters (2 from each session)
illusion_params <- list(
  list(col = "ca_amp",         label = "Calcium amplitude (a.u.)",     session = "calcium"),
  list(col = "ca_CTD90_ms",      label = "CTD90 (ms)",                   session = "calcium"),
  list(col = "contr_amp",        label = "Contraction amplitude (a.u.)", session = "contraction"),
  list(col = "Bodipy_intensity", label = "Bodipy intensity (a.u.)",      session = "contraction")
)

illusion_plots <- lapply(illusion_params, function(par) {
  well_src <- if (par$session == "calcium") ca_well_avgs else contr_well_avgs
  cell_src <- if (par$session == "calcium") calcium_gp else contraction_gp

  # Compute shared y-limits from the single-cell distribution (0.5–99.5th pct)
  # Well averages sit within this range by definition (means are less extreme)
  cell_vals <- cell_src[[par$col]]
  cell_vals <- cell_vals[is.finite(cell_vals)]
  y_lo  <- quantile(cell_vals, 0.005, na.rm = TRUE)
  y_hi  <- quantile(cell_vals, 0.995, na.rm = TRUE)
  y_range <- y_hi - y_lo
  shared_ylim <- c(y_lo - y_range * 0.02, y_hi + y_range * 0.08)

  n_above <- sum(cell_vals > y_hi, na.rm = TRUE)

  # Panel A: Well averages (traditional approach)
  p_well <- ggplot(well_src, aes(x = Condition, y = .data[[par$col]], fill = Condition)) +
    geom_point(aes(color = Condition), size = 2.5, alpha = 0.8,
               position = position_jitter(width = 0.15)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.4,
                 color = "grey20") +
    coord_cartesian(ylim = shared_ylim) +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = "Well averages", x = NULL, y = par$label) +
    theme(legend.position = "none", plot.title = element_text(size = 10))

  # Panel B: Single-cell data (our approach) — raincloud
  p_cell <- ggplot(cell_src, aes(x = Condition, y = .data[[par$col]],
                                 fill = Condition, color = Condition)) +
    geom_rain(alpha = 0.5, rain.side = "r",
              violin.args  = list(trim = TRUE, alpha = 0.45, color = NA),
              boxplot.args = list(outlier.shape = NA, width = 0.10,
                                  alpha = 0.6, color = "grey25"),
              point.args   = list(alpha = 0.07, size = 0.2, color = "grey35")) +
    coord_cartesian(ylim = shared_ylim) +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_fill_manual(values  = group_colors, drop = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = "Single-cell data", x = NULL, y = par$label) +
    theme(legend.position = "none", plot.title = element_text(size = 10))

  if (n_above > 0) {
    p_cell <- p_cell +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("\u2191 %d above", n_above),
               hjust = 1.05, vjust = 1.4, size = 2.5, color = "grey50",
               fontface = "italic")
  }

  list(p_well = p_well, p_cell = p_cell)
})

p_illusion <- ggpubr::ggarrange(
  plotlist = lapply(illusion_plots, function(x)
    ggpubr::ggarrange(x$p_well, x$p_cell, ncol = 2, widths = c(1, 1))),
  nrow = 4)
ggsave(file.path(fig_dir, paste0("fig_wellavg_vs_singlecell_", analysis_date, ".pdf")),
       p_illusion, width = 14, height = 16, device = cairo_pdf)
cat("  Saved: fig_wellavg_vs_singlecell\n")

# ---- 7b. Within-well heterogeneity: beeswarm by well ----------------------
# Show individual cells within wells for a key parameter
cat("  7c: Within-well heterogeneity plot...\n")

# Pick one parameter (calcium amplitude — most intuitive)
hetero_param <- "ca_amp"
hetero_label <- "Calcium amplitude (a.u.)"

# For readability, sample up to 6 wells per condition x cell line
set.seed(42)
wells_to_show <- calcium_gp %>%
  distinct(well_uid, Condition, Cell_Line) %>%
  group_by(Condition, Cell_Line) %>%
  slice_sample(n = 6, replace = FALSE) %>%
  ungroup()

hetero_data <- calcium_gp %>%
  filter(well_uid %in% wells_to_show$well_uid)

# Well averages for overlay
hetero_well_avg <- hetero_data %>%
  group_by(well_uid, Condition, Cell_Line) %>%
  summarise(avg = mean(.data[[hetero_param]], na.rm = TRUE), .groups = "drop")

p_hetero <- ggplot(hetero_data,
                   aes(x = well_uid, y = .data[[hetero_param]], color = Condition)) +
  geom_jitter(width = 0.2, alpha = 0.25, size = 0.3) +
  geom_point(data = hetero_well_avg, aes(y = avg),
             shape = 18, size = 3, color = "black") +
  facet_wrap(~ Cell_Line, scales = "free_x", nrow = 2, labeller = cell_line_fig_labeller) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  theme_project() +
  labs(title = "Single-cell heterogeneity within wells",
       subtitle = "Each dot = one cell. Black diamond = well average.",
       x = NULL, y = hetero_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom")
ggsave(file.path(fig_dir, paste0("fig_within_well_heterogeneity_", analysis_date, ".pdf")),
       p_hetero, width = 14, height = 8, device = cairo_pdf)
cat("  Saved: fig_within_well_heterogeneity\n")

# ---- 7d. PCA on single-cell multivariate data -----------------------------
# Since calcium and contraction are from separate sessions, run PCA per session
cat("  7d: PCA on single-cell parameters...\n")

run_session_pca <- function(data, pca_cols, load_labels, session_label,
                            signal_class_col, all_data = NULL) {
  pca_data <- data %>%
    select(cell_uid, well_uid, Condition, Cell_Line, all_of(pca_cols)) %>%
    drop_na()

  if (nrow(pca_data) < 50) {
    cat("    Skipped PCA for", session_label, "(too few complete cases)\n")
    return(invisible(NULL))
  }

  pca_mat <- pca_data %>% select(all_of(pca_cols)) %>% scale()
  pca_fit <- prcomp(pca_mat, center = FALSE, scale. = FALSE)
  var_explained <- summary(pca_fit)$importance[2, 1:2] * 100

  pca_scores <- as.data.frame(pca_fit$x[, 1:2])
  pca_scores$Condition <- pca_data$Condition
  pca_scores$Cell_Line <- pca_data$Cell_Line

  loadings <- as.data.frame(pca_fit$rotation[, 1:2])
  loadings$variable <- rownames(loadings)
  loadings$label <- load_labels[loadings$variable]
  arrow_scale <- 3

  p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(alpha = 0.15, size = 0.3) +
    stat_ellipse(level = 0.68, linewidth = 0.8) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "grey30", linewidth = 0.4) +
    geom_text(data = loadings,
              aes(x = PC1 * arrow_scale * 1.12, y = PC2 * arrow_scale * 1.12, label = label),
              inherit.aes = FALSE, size = 2.5, color = "grey20") +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = paste0("PCA: ", session_label, " (single-cell)"),
         subtitle = "Each dot = one cell. Ellipses = 68% of each group.",
         x = sprintf("PC1 (%.1f%% var)", var_explained[1]),
         y = sprintf("PC2 (%.1f%% var)", var_explained[2])) +
    theme(legend.position = "bottom")

  fname <- paste0("fig_pca_", tolower(gsub("[^A-Za-z]", "_", session_label)),
                  "_", analysis_date, ".pdf")
  ggsave(file.path(fig_dir, fname), p_pca, width = 12, height = 6, device = cairo_pdf)
  cat("  Saved:", fname, "\n")

  # Signal class PCA (good_paced vs irregular) — uses all QC-pass data
  if (!is.null(all_data) && !is.null(signal_class_col)) {
    pca_all <- all_data %>%
      select(cell_uid, Condition, Cell_Line,
             signal_class = all_of(signal_class_col),
             all_of(pca_cols)) %>%
      drop_na()

    if (nrow(pca_all) > 50) {
      pca_mat_all <- pca_all %>% select(all_of(pca_cols)) %>% scale()
      pca_fit_all <- prcomp(pca_mat_all, center = FALSE, scale. = FALSE)
      var_exp_all <- summary(pca_fit_all)$importance[2, 1:2] * 100

      pca_sc_all <- as.data.frame(pca_fit_all$x[, 1:2])
      pca_sc_all$signal_class <- pca_all$signal_class
      pca_sc_all$Cell_Line <- pca_all$Cell_Line
      pca_sc_all$Condition <- pca_all$Condition

      class_colors <- c(good_paced = "#2166AC", irregular = "#D95F02")
      p_class <- ggplot(pca_sc_all, aes(x = PC1, y = PC2, color = signal_class)) +
        geom_point(alpha = 0.15, size = 0.3) +
        stat_ellipse(level = 0.68, linewidth = 0.8) +
        facet_grid(Cell_Line ~ Condition, labeller = labeller(Cell_Line = cell_line_fig_labeller)) +
        scale_color_manual(values = class_colors, name = "Signal class") +
        theme_project() +
        labs(title = paste0("PCA: good-paced vs irregular (", session_label, ")"),
             subtitle = "Irregular cells occupy distinct regions of parameter space",
             x = sprintf("PC1 (%.1f%% var)", var_exp_all[1]),
             y = sprintf("PC2 (%.1f%% var)", var_exp_all[2])) +
        theme(legend.position = "bottom")

      fname2 <- paste0("fig_pca_signal_class_", tolower(gsub("[^A-Za-z]", "_", session_label)),
                        "_", analysis_date, ".pdf")
      ggsave(file.path(fig_dir, fname2), p_class, width = 14, height = 8, device = cairo_pdf)
      cat("  Saved:", fname2, "\n")
    }
  }
}

# Session 1: Calcium + Morphology (Cell.Area + Nucleus.Area valid here)
s1_pca_cols <- c("ca_amp", "ca_CTD50_ms", "ca_CTD90_ms", "ca_rise_ms", "ca_decay_ms",
                 "Cell.Area", "Nucleus.Area", "NC_ratio")
s1_pca_labels <- c(ca_amp = "Ca amp", ca_CTD50_ms = "CTD50", ca_CTD90_ms = "CTD90",
                   ca_rise_ms = "Ca rise", ca_decay_ms = "Ca decay",
                   Cell.Area = "Cell area", Nucleus.Area = "Nuc area", NC_ratio = "N/C")
run_session_pca(calcium_gp, s1_pca_cols, s1_pca_labels, "Calcium + Morphology",
                "ca_signal_class", calcium_all_classes)

# Session 2: Contraction + Nucleus.Area + Bodipy (Cell.Area unreliable here)
s2_pca_cols <- c("contr_amp", "contr_rise_ms", "contr_decay_ms",
                 "Nucleus.Area", "Bodipy_intensity")
s2_pca_labels <- c(contr_amp = "Contr amp", contr_rise_ms = "Contr rise",
                   contr_decay_ms = "Contr decay",
                   Nucleus.Area = "Nuc area (S2)", Bodipy_intensity = "Bodipy")
run_session_pca(contraction_gp, s2_pca_cols, s2_pca_labels, "Contraction + Bodipy",
                "contr_signal_class", contraction_all_classes)

# ---- 7e. Well-average vs mixed-model statistical comparison ----------------
cat("  7e: Well-average vs single-cell stats comparison...\n")

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
  # after clustering runs — see section 7f continued below.
}

# =============================================================================
# Section 7f: Subpopulation clustering analysis
# =============================================================================
# Cluster cells by functional profiles to ask: does lipid overload create
# distinct subpopulations? This is a key argument for single-cell resolution.
# Run per session (calcium+morphology and contraction+Bodipy separately).
cat("\n--- Section 7f: Subpopulation clustering analysis ---\n")

# run_subpop_analysis() and plot_cluster_traces() are defined in
# ../0_setup/cluster_helpers.R and sourced by run_pipeline.R at startup.

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

# --- Standardization confirmation for clustering inputs ---
# `scale()` (default center = TRUE, scale = TRUE) applies z-score
# standardization: each column has mean 0 and sd 1 before k-means. This
# mirrors cluster_helpers.R::run_subpop_analysis() exactly.
cat("\n  Confirming z-score standardization of clustering inputs...\n")
ca_mat_chk  <- calcium_gp %>%
  select(all_of(s1_clust_cols)) %>% drop_na() %>% as.matrix() %>% scale()
con_mat_chk <- contraction_gp %>%
  select(all_of(s2_clust_cols)) %>% drop_na() %>% as.matrix() %>% scale()

scaling_check <- bind_rows(
  tibble(session = "Calcium",     variable = colnames(ca_mat_chk),
         mean = round(colMeans(ca_mat_chk), 6),
         sd   = round(apply(ca_mat_chk, 2, sd), 6)),
  tibble(session = "Contraction", variable = colnames(con_mat_chk),
         mean = round(colMeans(con_mat_chk), 6),
         sd   = round(apply(con_mat_chk, 2, sd), 6))
)
cluster_variables_used <- bind_rows(
  tibble(session = "Calcium",     variable = s1_clust_cols,
         order = seq_along(s1_clust_cols)),
  tibble(session = "Contraction", variable = s2_clust_cols,
         order = seq_along(s2_clust_cols))
)
cat("\n  Post-scale() column means and SDs (should be 0 and 1):\n")
print(as.data.frame(scaling_check), row.names = FALSE)
write.csv(relabel_cell_line_table(scaling_check),
          file.path(table_dir, paste0("cluster_scaling_check_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(cluster_variables_used),
          file.path(table_dir, paste0("cluster_variables_used_",
                                      analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: cluster_scaling_check + cluster_variables_used\n")

subpop_calcium <- run_subpop_analysis(
  calcium_gp, s1_clust_cols, s1_clust_labels,
  "Calcium", fig_dir, analysis_date,
  k = k_calcium
)
if (!is.null(ca_results_all)) {
  p_ca_traces <- plot_cluster_traces(
    subpop_calcium, ca_cell_params_all, ca_results_all,
    session_label = "Calcium", y_label = "Calcium signal (baseline-corrected)",
    n_per_cluster = 8, fig_dir = fig_dir, analysis_date = analysis_date
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

  fig_clust_dir <- file.path(fig_dir, "figure_calcium_clusters")
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
  k = k_contraction
)
if (!is.null(contr_results_all)) {
  p_contr_traces <- plot_cluster_traces(
    subpop_contraction, contr_cell_params_all, contr_results_all,
    session_label = "Contraction + Bodipy",
    y_label = "Contraction (baseline-corrected)",
    n_per_cluster = 8, fig_dir = fig_dir, analysis_date = analysis_date
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

  fig_clust_dir <- file.path(fig_dir, "figure_contraction_clusters")
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

# ---- 7f continued: Two-panel power comparison figure ----
# Panel A: well-average vs single-cell mixed model p-values per parameter.
# Panel B: within-cluster condition effect per parameter — linear mixed model
#          (well as random intercept) on single-cell values within each cluster.
#          Falls back to KW on per-well means when the mixed model cannot be
#          estimated (too few wells or singular fit), flagged in the console.
if (nrow(wellavg_vs_mixed) > 0) {

  # Same palette as cluster composition / PCA plots in cluster_helpers.R
  clust_pal <- c("#2166AC", "#D95F02", "#1B9E77", "#E7298A", "#66A61E", "#7570B3")

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
        if (used_fallback)
          cat(sprintf("    [fallback KW] %s / Cluster %s\n", mp$label, k))

        tibble(
          parameter     = mp$label,
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

  # Shared y-axis factor (parameter rows only, no cluster proportion rows)
  param_levels <- rev(wellavg_vs_mixed$parameter)

  # Journal theme for both panels
  pwr_theme <- theme_bw() +
    theme(
      axis.text        = element_text(size = 7),
      axis.title       = element_text(size = 9),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 7),
      strip.text       = element_text(size = 9),
      axis.line        = element_line(linewidth = 0.5),
      axis.ticks       = element_line(linewidth = 0.5),
      panel.border     = element_rect(linewidth = 0.5),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      legend.position  = "bottom",
      legend.key.size  = unit(3.5, "mm")
    )

  # --- Panel A: well-average vs single-cell mixed model ---
  comp_models <- wellavg_vs_mixed %>%
    select(parameter, p_wellavg, p_mixed) %>%
    pivot_longer(c(p_wellavg, p_mixed), names_to = "m", values_to = "pvalue") %>%
    filter(!is.na(pvalue)) %>%
    transmute(
      parameter,
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
    annotate("text", x = -log10(0.05) + 0.25, y = 0.6,
             label = "p=0.05", size = 2.2, color = "grey50") +
    scale_color_manual(
      values = c("Well-average LM"          = "#999999",
                 "Single-cell mixed model"   = "#984EA3"),
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
  cluster_method_labels <- sort(unique(wc_rows$method))
  cluster_nums_used <- as.integer(gsub("Cluster (\\d+)", "\\1", cluster_method_labels))
  clust_color_vals <- setNames(clust_pal[cluster_nums_used], cluster_method_labels)
  clust_shape_vals <- setNames(rep(18L, length(cluster_method_labels)),
                               cluster_method_labels)

  wc_plot_data <- wc_rows %>%
    mutate(param_order = factor(parameter, levels = param_levels))

  p_B <- ggplot(wc_plot_data,
                aes(x = neg_log_p, y = param_order,
                    color = method, shape = method)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    annotate("text", x = -log10(0.05) + 0.25, y = 0.6,
             label = "p=0.05", size = 2.2, color = "grey50") +
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

  fig_pwr_dir <- file.path(fig_dir, "figure_power_comparison")
  dir.create(fig_pwr_dir, recursive = TRUE, showWarnings = FALSE)

  n_rows <- length(param_levels)
  fig_h_mm <- 25 + n_rows * 6.5

  ggsave(file.path(fig_pwr_dir, "figure_power_comparison.pdf"),
         fig_pwr,
         width = 183, height = fig_h_mm, units = "mm", device = "pdf")
  ggsave(file.path(fig_pwr_dir, "figure_power_comparison_A.pdf"),
         p_A + pwr_theme,
         width = 105, height = fig_h_mm, units = "mm", device = "pdf")
  ggsave(file.path(fig_pwr_dir, "figure_power_comparison_B.pdf"),
         p_B + pwr_theme + theme(axis.text.y = element_text(size = 7),
                                 axis.ticks.y = element_line(linewidth = 0.5)),
         width = 78, height = fig_h_mm, units = "mm", device = "pdf")

  cat(sprintf("  Saved: figure_power_comparison (A+B, %.0f mm tall)\n", fig_h_mm))
  if (any(wc_rows$used_fallback, na.rm = TRUE)) {
    cat("  Note: some clusters used KW fallback (see [fallback KW] lines above)\n")
  }
}

# =============================================================================
# Section 8: Export results (session-aware)
# =============================================================================
cat("\n--- Section 8: Exporting results ---\n")

# Good-paced cells per session
write.csv(relabel_cell_line_table(calcium_gp),
          file.path(table_dir, paste0("calcium_good_paced_", analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(contraction_gp),
          file.path(table_dir, paste0("contraction_good_paced_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: calcium_good_paced + contraction_good_paced CSVs\n")

# All QC-pass cells per session (includes irregular) — for auditing / QC review
write.csv(relabel_cell_line_table(calcium_all_classes),
          file.path(table_dir, paste0("calcium_all_qcpass_", analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(contraction_all_classes),
          file.path(table_dir, paste0("contraction_all_qcpass_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: calcium_all_qcpass + contraction_all_qcpass CSVs\n")

# Well-level merged (cross-session)
write.csv(relabel_cell_line_table(well_merged),
          file.path(table_dir, paste0("well_merged_cross_session_", analysis_date, ".csv")),
          row.names = FALSE)
cat("  Saved: well_merged_cross_session CSV\n")

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
if (exists("scaling_check"))        sheets[["Cluster_Scaling"]]    <- list(df = scaling_check,            p_cols = NULL)
if (exists("cluster_variables_used")) sheets[["Cluster_Variables"]] <- list(df = cluster_variables_used,  p_cols = NULL)
if (exists("well_counts_overall"))  sheets[["WellCounts_Summary"]] <- list(df = well_counts_overall,     p_cols = NULL)
if (exists("well_counts_detail"))   sheets[["WellCounts_Detail"]]  <- list(df = well_counts_detail,      p_cols = NULL)
if (exists("irregular_tests"))      sheets[["Irregular_Tests"]]    <- list(df = irregular_tests,         p_cols = c("LRT_p", "chisq_p"))
if (exists("irregular_glmm_fe"))    sheets[["Irregular_GLMM_FE"]]  <- list(df = irregular_glmm_fe,       p_cols = "p.value")
if (exists("ca_1624_ctrl_by_well")) sheets[["Irreg_1624_ctrl"]]    <- list(df = ca_1624_ctrl_by_well,    p_cols = NULL)
if (exists("well_cross_cor"))       sheets[["Well_Cross_Cor"]]     <- list(df = well_cross_cor,          p_cols = "p_value")

# Well-level QC report (from run_pipeline.R well_qc_gate tables)
if (exists("ca_well_qc_gate") && exists("contr_well_qc_gate")) {
  well_qc_report <- bind_rows(ca_well_qc_gate, contr_well_qc_gate) %>%
    arrange(modality, well_qc_pass, well_uid) %>%
    select(modality, well_uid, Plate, Well, Condition, Cell_Line,
           n_total, n_good_paced, n_irregular, pct_good_paced,
           pct_irregular, well_qc_pass, qc_fail_reason)
  sheets[["Well_QC"]] <- list(df = well_qc_report, p_cols = NULL)
}

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
cat("Output files:\n")
cat("  Figures directory:", fig_dir, "\n")
cat("    -", length(violin_specs), "raincloud plots with p-values and n per group (good_paced cells)\n")
cat("    - 6 correlation heatmaps (within-session + cross-session)\n")
cat("    - 8 scatter plots (within-session single-cell + cross-session well-level)\n")
cat("    - 1 irregular cell proportion barplot (per session)\n")
cat("    Single-cell manuscript figures:\n")
cat("    - fig_wellavg_vs_singlecell (well avg illusion)\n")
cat("    - fig_variance_decomposition (ICC stacked bar)\n")
cat("    - fig_within_well_heterogeneity (beeswarm by well)\n")
cat("    - fig_pca_calcium_morphology + fig_pca_contraction_bodipy\n")
cat("    - fig_pca_signal_class (per session)\n")
cat("    - fig_power_comparison (well-avg vs mixed model)\n")
cat("    - fig_cluster_profiles + fig_pca_clusters (subpopulation analysis)\n")
cat("    - fig_cluster_composition (condition-specific subpopulation shifts)\n")
cat("  Tables directory:", table_dir, "\n")
cat("    - calcium_good_paced + contraction_good_paced (session-specific)\n")
cat("    - calcium_all_qcpass + contraction_all_qcpass (incl. irregular)\n")
cat("    - irregular_cells_calcium + irregular_cells_contraction\n")
cat("    - well_merged_cross_session (well-level averages)\n")
cat("    - summary_statistics (per session)\n")
cat("    - statistical_results.xlsx (Results_Summary | ANOVA | Posthoc |\n")
cat("        Fixed_Effects | Summary_Calcium | Summary_Contraction |\n")
cat("        Irregular_Summary | Variance_Decomp | Power_Comparison)\n")
cat("================================================\n")



# =============================================================================
# Section 9: Figure 2 — variance decomposition + well-avg vs single-cell
#
#   A  Variance decomposition (ICC) — all measured parameters
#   B  Ca amplitude, CTD50, Contr. rise, Contr. decay (well avg | single-cell)
#
# Output: 3_output/Figures/figure2/
# (Re-loads calcium_gp / contraction_gp from the CSVs just written in Section 3
#  so factor levels are explicit. OK here because no downstream section uses them.)
# =============================================================================

fig_out <- file.path(main_dir, "Figures", "figure2")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

# ---- Load processed data ----------------------------------------------------
load_latest <- function(dir, stem) {
  files <- list.files(dir, pattern = paste0("^", stem), full.names = TRUE)
  if (length(files) == 0) stop("No file matching '", stem, "' in ", dir)
  files[order(file.mtime(files), decreasing = TRUE)][1]
}

calcium_gp     <- read_csv(load_latest(table_dir, "calcium_good_paced"),
                            show_col_types = FALSE)
contraction_gp <- read_csv(load_latest(table_dir, "contraction_good_paced"),
                            show_col_types = FALSE)

set_factors <- function(df) {
  df$Condition <- factor(df$Condition, levels = group_levels)
  df$Cell_Line <- factor(df$Cell_Line, levels = c("1.6", "1.6.24"))
  df
}
calcium_gp     <- set_factors(calcium_gp)
contraction_gp <- set_factors(contraction_gp)

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
# Exact same geom_rain style as fig_wellavg_vs_singlecell in 4_plot_and_stats.R
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
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
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
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
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

if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)

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

fig2 <- plot_grid(p_A, panel_B,
                  nrow        = 2,
                  rel_heights = c(1, 2.6),
                  labels      = c("A", "B"),
                  label_size  = 9,
                  label_fontface = "bold")

# =============================================================================
# Save
# =============================================================================
cat("Saving to:", fig_out, "\n")

ggsave(file.path(fig_out, "figure2.pdf"), fig2,
       width = 183, height = 247, units = "mm", device = "pdf")

ggsave(file.path(fig_out, "figure2_A.pdf"),
       p_A + theme(legend.position = "bottom"),
       width = 183, height = 80, units = "mm", device = "pdf")

ggsave(file.path(fig_out, "figure2_B.pdf"),
       panel_B,
       width = 183, height = 160, units = "mm", device = "pdf")

cat("Done. Saved:\n")
cat("  figure2.pdf (183 x 247 mm)\n")
cat("  figure2_A.pdf, figure2_B.pdf\n")
cat("\nREMINDER for figure caption:\n")
cat("  Well-average panels: points = individual wells, crossbar = mean\n")
cat("  Single-cell panels: violin = KDE, box = median + IQR, dots = cells\n")

# =============================================================================
# Append Variance_Decomp sheet to statistical_results.xlsx
# =============================================================================
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

xlsx_files <- list.files(table_dir,
                         pattern = "^statistical_results.*\\.xlsx$",
                         full.names = TRUE)
if (length(xlsx_files) > 0) {
  xlsx_path <- xlsx_files[order(file.mtime(xlsx_files), decreasing = TRUE)][1]
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
  cat("  WARNING: No statistical_results.xlsx found — Variance_Decomp not saved.\n")
}



# =============================================================================
# Section 10: Graphical-abstract trace illustration
#
# 8 individual single-cell traces per modality (calcium + contraction),
# small-multiples layout (4 x 2 grid). No axis labels.
# Output: 3_output/Figures/figure_trace_illustration/
# =============================================================================
main_dir      <- "../3_output"
analysis_date <- format(Sys.Date(), "%Y-%m-%d")
n_cells       <- 8
n_col         <- 4          # columns in facet grid  → 4 x 2
trace_window  <- c(0, 5000) # ms

# ---- Load ----
cat("Loading RDS files...\n")
ca_params   <- readRDS(file.path(main_dir, "processed", "calcium_all_cells.rds"))
ca_results  <- readRDS(file.path(main_dir, "processed", "calcium_results_all.rds"))
con_params  <- readRDS(file.path(main_dir, "processed", "contraction_all_cells.rds"))
con_results <- readRDS(file.path(main_dir, "processed", "contraction_results_all.rds"))

# Pick a random well with enough good-paced cells in both modalities
candidate_wells <- intersect(
  ca_params$well_uid[ca_params$qc_pass & ca_params$signal_class == "good_paced"],
  con_params$well_uid[con_params$qc_pass & con_params$signal_class == "good_paced"]
)
well_counts <- table(candidate_wells)
eligible    <- names(well_counts[well_counts >= n_cells])
if (length(eligible) == 0) eligible <- names(well_counts)  # relax if needed
target_well <- sample(eligible, 1)
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
trace_fig_dir <- file.path(main_dir, "Figures", "figure_trace_illustration")
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
