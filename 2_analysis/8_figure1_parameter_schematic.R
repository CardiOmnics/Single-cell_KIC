# =============================================================================
# 8_figure1_parameter_schematic.R
# Schematic of the functional parameters, drawn on one real Ca2+ transient and
# one real contraction beat (proposed for Fig. 1 or a supplementary figure).
#
#   A  Ca2+ transient: peak amplitude, rise time, decay time, CTD50, CTD90
#   B  Contraction:    amplitude, rise time, decay time
#
# Each panel shows one beat of one control cell whose parameters lie closest
# to the control median (z-scored distance), so the example is typical rather
# than hand-picked. Within that cell, the beat shown is the one whose per-beat
# values lie closest to the cell's own medians (the medians are what enter the
# analysis). The marked intervals are computed with the same functions that
# produce the analysed parameters (0_setup/contraction_trace_utils.R):
#   Ca rise time     last sample <= 10% of peak before the peak -> peak
#   Ca decay time    peak -> first sample <= 20% of peak (80% recovery)
#   CTD50 / CTD90    upstroke -> downstroke crossing at 50% / 10% of peak
#   Contr. rise time last sample <= 20% of peak before the peak -> peak
#   Contr. decay     peak -> first sample <= 20% of peak (80% relaxation)
# The y axis starts at the zero of the analysed trace. For calcium this is the
# CyteSeer trimmed trace without further baseline correction, so amplitude and
# the 10/20/50% levels are fractions of the peak value above that zero,
# exactly as in the analysis.
#
# Requires (from run_pipeline.R): ca_cell_params_all, ca_results_all,
#   contr_cell_params_all, contr_results_all (well QC applied to qc_pass).
#   Standalone: loads the cached RDS files and keeps only the cells listed in
#   the newest *_good_paced table.
# Output: 3_output/Figures/Publication/figure_parameter_schematic/
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

if (!exists("main_dir"))      main_dir      <- "../3_output"
if (!exists("analysis_date")) analysis_date <- format(Sys.Date(), "%Y-%m-%d")
if (!exists("pub_dir"))       pub_dir       <- file.path(main_dir, "Figures", "Publication")
if (!exists("calc_ctd"))      source("../0_setup/contraction_trace_utils.R")
if (!exists("journal_theme")) source("../0_setup/theme.R")
if (!exists("relabel_condition")) source("../0_setup/colors.R")

# ---- Data -------------------------------------------------------------------
if (!exists("ca_results_all") || is.null(ca_results_all) ||
    !exists("contr_results_all") || is.null(contr_results_all)) {
  cache_dir <- if (exists("processed_dir")) processed_dir else "../3_output/processed"
  ca_cell_params_all    <- relabel_condition(readRDS(file.path(cache_dir, "calcium_all_cells.rds")))
  ca_results_all        <- readRDS(file.path(cache_dir, "calcium_results_all.rds"))
  contr_cell_params_all <- relabel_condition(readRDS(file.path(cache_dir, "contraction_all_cells.rds")))
  contr_results_all     <- readRDS(file.path(cache_dir, "contraction_results_all.rds"))
  newest <- function(stem) {
    f <- list.files(file.path(main_dir, "Tables"), pattern = paste0("^", stem, "_[0-9-]+\\.csv$"),
                    full.names = TRUE)
    f[which.max(file.mtime(f))]
  }
  keep_ca <- read_csv(newest("calcium_good_paced"),     show_col_types = FALSE)$cell_uid
  keep_co <- read_csv(newest("contraction_good_paced"), show_col_types = FALSE)$cell_uid
  ca_cell_params_all$qc_pass    <- ca_cell_params_all$cell_uid    %in% keep_ca
  contr_cell_params_all$qc_pass <- contr_cell_params_all$cell_uid %in% keep_co
}

# Control cell closest to the control median on the given parameters
pick_typical_cell <- function(params, cols, min_beats) {
  ctrl <- params %>%
    mutate(row_idx = row_number()) %>%
    filter(qc_pass, signal_class == "good_paced", Condition == "Control",
           n_beats >= min_beats) %>%
    drop_na(all_of(cols))
  z <- scale(as.matrix(ctrl[, cols]))
  med <- apply(z, 2, median)
  d <- sqrt(rowSums(sweep(z, 2, med)^2))
  ctrl[which.min(d), ]
}

# One beat of the trace with its landmarks, measured as in the pipeline
pick_beat <- function(res, rise_frac, decay_frac, ctd = FALSE, pad_ms = 150) {
  df <- as.data.frame(res$trace); pk <- res$peaks
  df <- df[!is.na(df$signal_bc) & !is.na(df$time_ms), ]
  kin <- lapply(seq_len(nrow(pk)), function(i)
    extract_beat_kinetics(df, pk[i, ], rise_frac = rise_frac, decay_frac = decay_frac))
  beats <- tibble(rise  = map_dbl(kin, ~ if (is.null(.x)) NA_real_ else .x$rise_time_ms),
                  decay = map_dbl(kin, ~ if (is.null(.x)) NA_real_ else .x$decay_time_ms))
  if (ctd) beats <- beats %>% mutate(ctd50 = calc_ctd(df, pk, 0.5), ctd90 = calc_ctd(df, pk, 0.9))
  meds  <- map_dbl(beats, median, na.rm = TRUE)
  score <- rowSums(abs(sweep(as.matrix(beats), 2, meds)) / meds)
  b <- which.min(score)

  amp <- pk$peak_amp[b]; t_pk <- pk$peak_time[b]
  lm <- list(amp = amp, t_peak = 0,
             t_rise0 = -beats$rise[b], t_decay1 = beats$decay[b],
             rise_ms = beats$rise[b], decay_ms = beats$decay[b])
  if (ctd) {
    cross <- function(frac) {          # same crossing rule as calc_ctd()
      thr <- (1 - frac) * amp; pi <- pk$peak_idx[b]
      up <- tail(which(df$signal_bc[1:pi] <= thr), 1)
      dn <- pi + which(df$signal_bc[pi:nrow(df)] <= thr)[1] - 1
      c(df$time_ms[up], df$time_ms[dn]) - t_pk
    }
    lm$ctd50 <- cross(0.5); lm$ctd90 <- cross(0.9)
    stopifnot(isTRUE(all.equal(diff(lm$ctd50), beats$ctd50[b])),
              isTRUE(all.equal(diff(lm$ctd90), beats$ctd90[b])))
  }
  t0 <- min(c(lm$t_rise0, lm$ctd90[1])) - pad_ms
  t1 <- max(c(lm$t_decay1, lm$ctd90[2])) + pad_ms
  win <- df %>% transmute(time_ms = time_ms - t_pk, signal_bc) %>%
    filter(time_ms >= t0, time_ms <= t1)
  list(trace = win, lm = lm, beat = b, n_beats = nrow(pk), meds = meds)
}

# ---- Pick cells and beats ---------------------------------------------------
ca_cell <- pick_typical_cell(ca_cell_params_all,
                             c("mean_amp", "CTD50_ms", "CTD90_ms", "rise_time_ms", "decay_time_ms"),
                             min_beats = 8)
co_cell <- pick_typical_cell(contr_cell_params_all,
                             c("mean_amp", "rise_time_ms", "decay_time_ms"),
                             min_beats = 4)
cat(sprintf("  Ca example:          %s (%s)\n", ca_cell$cell_uid, ca_cell$Cell_Line))
cat(sprintf("  Contraction example: %s (%s)\n", co_cell$cell_uid, co_cell$Cell_Line))

ca_pick <- pick_beat(ca_results_all[[ca_cell$row_idx]], rise_frac = 0.1, decay_frac = 0.2, ctd = TRUE)
co_pick <- pick_beat(contr_results_all[[co_cell$row_idx]], rise_frac = 0.2, decay_frac = 0.2)
ca_beat <- ca_pick$trace; ca_lm <- ca_pick$lm
co_beat <- co_pick$trace; co_lm <- co_pick$lm

# ---- Drawing helpers ----------------------------------------------------------
ink   <- "grey20"
guide <- "grey60"
arr   <- arrow(length = unit(1.2, "mm"), ends = "both", type = "closed")
txt   <- 2.3   # mm, ~6.5 pt

h_interval <- function(x0, x1, y, label, vjust = -0.5) {
  list(annotate("segment", x = x0, xend = x1, y = y, yend = y,
                colour = ink, linewidth = 0.35, arrow = arr),
       annotate("text", x = (x0 + x1) / 2, y = y, label = label,
                vjust = vjust, size = txt, colour = ink, lineheight = 0.9))
}
level_line <- function(x0, x1, y) {
  annotate("segment", x = x0, xend = x1, y = y, yend = y,
           colour = guide, linewidth = 0.25, linetype = "dashed")
}
v_tick <- function(x, y0, y1) {
  annotate("segment", x = x, xend = x, y = y0, yend = y1,
           colour = guide, linewidth = 0.25, linetype = "dotted")
}

base_plot <- function(beat, colour, y_lab) {
  ggplot(beat, aes(time_ms, signal_bc)) +
    geom_hline(yintercept = 0, colour = guide, linewidth = 0.25) +
    geom_line(colour = colour, linewidth = 0.8) +
    geom_point(colour = colour, size = 0.5) +
    labs(x = "Time from peak (ms)", y = y_lab) +
    journal_theme +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 9, face = "bold"))
}

# ---- Panel A: calcium ---------------------------------------------------------
a  <- ca_lm$amp
x_lo <- min(ca_beat$time_ms); x_hi <- max(ca_beat$time_ms)
p_ca <- base_plot(ca_beat, "#1B6CA8", expression(Ca^{2*"+"}~signal~"(a.u.)")) +
  # amplitude
  annotate("segment", x = x_lo + 45, xend = x_lo + 45, y = 0, yend = a,
           colour = ink, linewidth = 0.35, arrow = arr) +
  annotate("text", x = x_lo + 22, y = a / 2, label = "Peak amplitude", angle = 90,
           vjust = 0, size = txt, colour = ink) +
  level_line(x_lo + 20, 0, a) +
  # CTD50 and CTD90 (width at 50% and 10% of peak)
  level_line(x_lo, x_hi, 0.5 * a) + level_line(x_lo, x_hi, 0.1 * a) +
  h_interval(ca_lm$ctd50[1], ca_lm$ctd50[2], 0.5 * a, "CTD50") +
  h_interval(ca_lm$ctd90[1], ca_lm$ctd90[2], 0.1 * a, "CTD90") +
  annotate("text", x = x_hi, y = c(0.5, 0.1) * a, label = c("50%", "10%"),
           hjust = 1, vjust = -0.4, size = 2, colour = guide) +
  # rise (10% -> peak) and decay (peak -> 20%) drawn above the trace
  v_tick(ca_lm$t_rise0, 0.1 * a, 1.16 * a) + v_tick(0, a, 1.40 * a) +
  v_tick(ca_lm$t_decay1, 0.2 * a, 1.40 * a) +
  h_interval(ca_lm$t_rise0, 0, 1.16 * a, "Rise time\n(10% to peak)") +
  h_interval(0, ca_lm$t_decay1, 1.40 * a, "Decay time\n(peak to 80% recovery)") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.02))) +
  coord_cartesian(ylim = c(0, 1.65 * a), clip = "off") +
  labs(title = expression(bold(Ca^{2*"+"}~transient)))

# ---- Panel B: contraction ----------------------------------------------------
b <- co_lm$amp
xb_lo <- min(co_beat$time_ms); xb_hi <- max(co_beat$time_ms)
p_co <- base_plot(co_beat, "#C94040", "Contraction signal (a.u.)") +
  annotate("segment", x = xb_lo + 35, xend = xb_lo + 35, y = 0, yend = b,
           colour = ink, linewidth = 0.35, arrow = arr) +
  annotate("text", x = xb_lo + 15, y = b / 2, label = "Amplitude", angle = 90,
           vjust = 0, size = txt, colour = ink) +
  level_line(xb_lo + 15, 0, b) +
  level_line(xb_lo, xb_hi, 0.2 * b) +
  annotate("text", x = xb_hi, y = 0.2 * b, label = "20%", hjust = 1, vjust = -0.4,
           size = 2, colour = guide) +
  v_tick(co_lm$t_rise0, 0.2 * b, 1.16 * b) + v_tick(0, b, 1.40 * b) +
  v_tick(co_lm$t_decay1, 0.2 * b, 1.40 * b) +
  h_interval(co_lm$t_rise0, 0, 1.16 * b, "Rise time\n(20% to peak)") +
  h_interval(0, co_lm$t_decay1, 1.40 * b, "Decay time\n(peak to 80% relaxation)") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.02))) +
  coord_cartesian(ylim = c(0, 1.65 * b), clip = "off") +
  labs(title = "Contraction")

if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
fig <- cowplot::plot_grid(p_ca, p_co, ncol = 2, rel_widths = c(1.15, 1),
                          labels = c("A", "B"), label_size = 9)

out_dir <- file.path(pub_dir, "figure_parameter_schematic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, paste0("figure_parameter_schematic_", analysis_date, ".pdf")),
       fig, width = 183, height = 75, units = "mm", device = "pdf")
ggsave(file.path(out_dir, paste0("figure_parameter_schematic_", analysis_date, ".png")),
       fig, width = 183, height = 75, units = "mm", dpi = 300)

write_csv(tibble(
  panel      = c("Ca", "Contraction"),
  cell_uid   = c(ca_cell$cell_uid, co_cell$cell_uid),
  beat_shown = c(ca_pick$beat, co_pick$beat),
  n_beats    = c(ca_pick$n_beats, co_pick$n_beats),
  amplitude  = c(ca_lm$amp, co_lm$amp),
  rise_ms    = c(ca_lm$rise_ms, co_lm$rise_ms),
  decay_ms   = c(ca_lm$decay_ms, co_lm$decay_ms),
  CTD50_ms   = c(diff(ca_lm$ctd50), NA),
  CTD90_ms   = c(diff(ca_lm$ctd90), NA)),
  file.path(out_dir, paste0("figure_parameter_schematic_cells_", analysis_date, ".csv")))
cat("  Saved:", out_dir, "\n")
