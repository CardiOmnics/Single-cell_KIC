# ============================================================
# Master pipeline: Single-cell calcium, contraction & Bodipy
# Run from top to bottom — no interactive prompts.
# To exclude specific wells, edit the vectors below before running.
# ============================================================

rm(list = ls())
set.seed(123)

# Use cloud CRAN mirror so install.packages() never prompts for selection
options(repos = c(CRAN = "https://cloud.r-project.org"))

core_packages <- c("tidyverse", "ggpubr", "broom", "lme4", "emmeans", "openxlsx",
                   "lmerTest", "broom.mixed", "pracma")
missing_packages <- core_packages[!vapply(core_packages, requireNamespace,
                                          logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) install.packages(missing_packages)
suppressPackageStartupMessages(
  invisible(lapply(core_packages, library, character.only = TRUE))
)

# Run this script with 2_analysis/ as the working directory, e.g.
#   setwd("<path to repository>/2_analysis"); source("run_pipeline.R")
# All paths below are relative to 2_analysis/.

# Global settings
#   main_dir     = 3_output/                     all pipeline outputs
#   processed    = 3_output/processed/           per-cell RDS files (Zenodo)
#   pub_dir      = 3_output/Figures/Publication  main-figure panels
#   expl_dir     = 3_output/Figures/Exploratory  supplementary figures,
#                  single panels, diagnostics, trace galleries
main_dir      <- "../3_output"
analysis_date <- format(Sys.Date(), "%Y-%m-%d")

processed_dir <- file.path(main_dir, "processed")
pub_dir       <- file.path(main_dir, "Figures", "Publication")
expl_dir      <- file.path(main_dir, "Figures", "Exploratory")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pub_dir,       showWarnings = FALSE, recursive = TRUE)
dir.create(expl_dir,      showWarnings = FALSE, recursive = TRUE)

# Load project-wide settings
source("../0_setup/colors.R")
source("../0_setup/theme.R")
source("../0_setup/helpers.R")
source("../0_setup/contraction_trace_utils.R")

theme_set(theme_project())

cat("========================================\n")
cat("  Single-cell analysis pipeline\n")
cat("  Date:", analysis_date, "\n")
cat("========================================\n\n")

# ============================================================
# Steps to run — set to FALSE to skip
# ============================================================
# run_contraction / run_calcium = TRUE re-process the raw CyteSeer exports
# (download from Zenodo first, see 1_data/README.md) and write the RDS files
# in 3_output/processed/. FALSE uses the RDS files from Zenodo placed in
# 3_output/processed/, which is all that is needed to reproduce every figure,
# table and statistic.
run_contraction <- FALSE
run_calcium     <- FALSE
run_qc_plots    <- FALSE  # QC diagnostics + trace galleries (not in the paper)
run_figures     <- TRUE   # all manuscript figures, tables and statistics

# ============================================================
# Well exclusions — edit these before running if needed
# Format: "plate_well", e.g. "1.6_B02" or "1.6.24_D02"
# ============================================================
contr_exclude_wells <- c(
  "1.6.24_B04",  # beating at ~2 Hz (double pacing freq); no cells follow 1 Hz pacing
  "1.6_B03"      # beating at ~2 Hz (double pacing freq); no cells follow 1 Hz pacing
)
# Manual exclusions are reserved for PACING failures, which no automatic gate
# detects: these wells beat at ~2 Hz (double the 1 Hz pacing frequency), so no
# cell follows the pacemaker. Everything else is handled by the automatic
# well-level gates below (min_good_paced + responder rate).
# (1.6.24_E04, a segmentation artifact, is caught by the responder-rate gate.)
ca_exclude_wells    <- c(
  "1.6.24_E06",  # beating at ~2 Hz (double pacing freq); no cells follow 1 Hz pacing
  "1.6_C05"      # beating at ~2 Hz (double pacing freq); no cells follow 1 Hz pacing
)

# ============================================================
# Well-level QC: minimum number of good_paced cells per well
# Wells below this threshold are excluded from main analyses.
# ============================================================
min_good_paced <- 10


# ============================================================
# STEP 1: Process contraction data + Bodipy + Nucleus.Area
# ============================================================

if (run_contraction) {
  cat("STEP 1: Processing contraction data...\n")
  cat("────────────────────────────────────────\n")
  source("1_process_contraction.R")

  cat("\nContraction QC summary (sorted by % passing):\n\n")
  print(as.data.frame(contr_well_qc_summary), row.names = FALSE)
} else {
  cat("STEP 1: Skipped — loading saved data from disk\n")
  rds <- file.path(processed_dir, "contraction_all_cells.rds")
  if (!file.exists(rds)) stop("No contraction RDS in ", processed_dir, ". Download it from Zenodo ",
                              "(see 1_data/README.md) or set run_contraction = TRUE.")
  contr_cell_params_all <<- readRDS(rds)
  cat("  Loaded:", rds, "\n")
  results_rds <- file.path(processed_dir, "contraction_results_all.rds")
  if (file.exists(results_rds)) {
    contr_results_all <<- readRDS(results_rds)
    cat("  Loaded:", results_rds, "\n\n")
  } else {
    contr_results_all <<- NULL   # not yet saved; cluster traces skipped
    cat("  (contraction_results_all.rds not found —",
        "run step 1 once to save it)\n\n")
  }
}

# Apply contraction well exclusions (runs whether data was freshly
# processed or loaded from disk, so exclusions are always enforced)
if (length(contr_exclude_wells) > 0) {
  n_before <- sum(contr_cell_params_all$qc_pass)
  excl <- contr_cell_params_all$well_uid %in% contr_exclude_wells
  contr_cell_params_all$qc_pass[excl] <- FALSE
  n_after <- sum(contr_cell_params_all$qc_pass)
  cat(sprintf("  -> Excluded %d contraction wells, removed %d additional cells\n",
              length(contr_exclude_wells), n_before - n_after))
}

# Well-level QC: exclude wells with fewer than min_good_paced good_paced cells
contr_well_qc_gate <- contr_cell_params_all %>%
  group_by(well_uid, Plate, Well, Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = sum(signal_class == "irregular"),
    n_manually_excluded = sum(!qc_pass),
    pct_good_paced = round(100 * n_good_paced / n_total, 1),
    pct_irregular  = round(100 * n_irregular / n_total, 1),
    .groups = "drop"
  ) %>%
  mutate(
    well_qc_pass = n_good_paced >= min_good_paced,
    qc_fail_reason = case_when(
      !well_qc_pass ~ sprintf("n_good_paced = %d (< %d)",
                               n_good_paced, min_good_paced),
      n_manually_excluded == n_total ~ "manually excluded",
      TRUE ~ NA_character_
    ),
    modality = "Contraction"
  )
contr_failed_wells <- contr_well_qc_gate %>% filter(!well_qc_pass)
if (nrow(contr_failed_wells) > 0) {
  n_before <- sum(contr_cell_params_all$qc_pass)
  contr_cell_params_all$qc_pass[
    contr_cell_params_all$well_uid %in% contr_failed_wells$well_uid] <- FALSE
  n_after <- sum(contr_cell_params_all$qc_pass)
  cat(sprintf("  -> Well-level QC: excluded %d contraction wells with < %d good_paced cells, removed %d cells\n",
              nrow(contr_failed_wells), min_good_paced, n_before - n_after))
  cat("     Failed wells:", paste(contr_failed_wells$well_uid, collapse = ", "), "\n")
}

n_pass  <- sum(contr_cell_params_all$qc_pass)
n_total <- nrow(contr_cell_params_all)
n_good  <- sum(contr_cell_params_all$qc_pass & contr_cell_params_all$signal_class == "good_paced")
n_irreg <- sum(contr_cell_params_all$qc_pass & contr_cell_params_all$signal_class == "irregular")
cat(sprintf("  -> Contraction: %d / %d cells passing QC (%.1f%%)\n",
            n_pass, n_total, 100 * n_pass / n_total))
cat(sprintf("     good_paced: %d | irregular: %d\n\n", n_good, n_irreg))


# ============================================================
# STEP 2: Process calcium data + Cell.Area + Nucleus.Area
# ============================================================

if (run_calcium) {
  cat("STEP 2: Processing calcium data...\n")
  cat("────────────────────────────────────────\n")
  source("2_process_calcium.R")

  cat("\nCalcium QC summary (sorted by % passing):\n\n")
  print(as.data.frame(ca_well_qc_summary), row.names = FALSE)
} else {
  cat("STEP 2: Skipped — loading saved data from disk\n")
  rds <- file.path(processed_dir, "calcium_all_cells.rds")
  if (!file.exists(rds)) stop("No calcium RDS in ", processed_dir, ". Download it from Zenodo ",
                              "(see 1_data/README.md) or set run_calcium = TRUE.")
  ca_cell_params_all <<- readRDS(rds)
  cat("  Loaded:", rds, "\n")
  results_rds <- file.path(processed_dir, "calcium_results_all.rds")
  if (file.exists(results_rds)) {
    ca_results_all <<- readRDS(results_rds)
    cat("  Loaded:", results_rds, "\n\n")
  } else {
    ca_results_all <<- NULL   # not yet saved; cluster traces skipped
    cat("  (calcium_results_all.rds not found —",
        "run step 2 once to save it)\n\n")
  }
}

# Apply calcium well exclusions (runs whether data was freshly
# processed or loaded from disk, so exclusions are always enforced)
if (length(ca_exclude_wells) > 0) {
  n_before <- sum(ca_cell_params_all$qc_pass)
  excl <- ca_cell_params_all$well_uid %in% ca_exclude_wells
  ca_cell_params_all$qc_pass[excl] <- FALSE
  n_after <- sum(ca_cell_params_all$qc_pass)
  cat(sprintf("  -> Excluded %d calcium wells, removed %d additional cells\n",
              length(ca_exclude_wells), n_before - n_after))
}

# ---- Well-level signal-quality gate: responder rate -------------------------
# Fraction of segmented cells that CyteSeer flags as having a detectable
# calcium transient ("Ci on Cm API is Responder"). This is computed by the
# instrument software, independently of the trace analysis in
# 2_process_calcium.R, so it is a measure of acquisition/signal quality rather
# than of the biological readouts being modelled.
#
# It is needed because the cell-level AC1 gate scores RHYTHMICITY, not signal
# strength: a weak but periodic trace has high autocorrelation and passes,
# so a well can be ~100% "good_paced" while most of its cells have no
# detectable transient at all.
#
# The cutoff is Otsu's method on the well-level distribution — the same
# data-driven estimator already used for the cell-level AC1 threshold, so there
# is no hand-chosen number and no well is named in the code.
#
# Contraction is not gated this way: that pipeline (LicensedPipelineAlgorithm)
# emits no equivalent responder flag.
ca_responder_by_well <- ca_cell_params_all %>%
  group_by(well_uid) %>%
  summarise(pct_responder = 100 * mean(is_responder, na.rm = TRUE),
            .groups = "drop")

resp_otsu <- otsu_threshold(ca_responder_by_well$pct_responder)
cat(sprintf("  -> Responder-rate Otsu threshold: %.1f%% (n = %d wells)\n",
            resp_otsu, nrow(ca_responder_by_well)))

# Well-level QC: minimum good_paced cells AND minimum responder rate
ca_well_qc_gate <- ca_cell_params_all %>%
  group_by(well_uid, Plate, Well, Condition, Cell_Line) %>%
  summarise(
    n_total      = n(),
    n_good_paced = sum(signal_class == "good_paced"),
    n_irregular  = sum(signal_class == "irregular"),
    n_manually_excluded = sum(!qc_pass),
    pct_good_paced = round(100 * n_good_paced / n_total, 1),
    pct_irregular  = round(100 * n_irregular / n_total, 1),
    pct_responder  = round(100 * mean(is_responder, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  mutate(
    pass_n_good   = n_good_paced >= min_good_paced,
    pass_responder = pct_responder >= resp_otsu,
    well_qc_pass  = pass_n_good & pass_responder,
    qc_fail_reason = case_when(
      !pass_n_good & !pass_responder ~ sprintf(
        "n_good_paced = %d (< %d); responders = %.1f%% (< %.1f%%)",
        n_good_paced, min_good_paced, pct_responder, resp_otsu),
      !pass_n_good ~ sprintf("n_good_paced = %d (< %d)",
                             n_good_paced, min_good_paced),
      !pass_responder ~ sprintf("responders = %.1f%% (< %.1f%% Otsu)",
                                pct_responder, resp_otsu),
      n_manually_excluded == n_total ~ "manually excluded",
      TRUE ~ NA_character_
    ),
    modality = "Calcium"
  )
ca_failed_wells <- ca_well_qc_gate %>% filter(!well_qc_pass)
if (nrow(ca_failed_wells) > 0) {
  n_before <- sum(ca_cell_params_all$qc_pass)
  ca_cell_params_all$qc_pass[
    ca_cell_params_all$well_uid %in% ca_failed_wells$well_uid] <- FALSE
  n_after <- sum(ca_cell_params_all$qc_pass)
  cat(sprintf("  -> Well-level QC: excluded %d calcium wells, removed %d cells\n",
              nrow(ca_failed_wells), n_before - n_after))
  cat(sprintf("     by n_good_paced < %d : %d well(s)\n",
              min_good_paced, sum(!ca_failed_wells$pass_n_good)))
  cat(sprintf("     by responders < %.1f%%: %d well(s)\n",
              resp_otsu, sum(!ca_failed_wells$pass_responder)))
  for (i in seq_len(nrow(ca_failed_wells))) {
    cat(sprintf("       %-12s %s\n", ca_failed_wells$well_uid[i],
                ca_failed_wells$qc_fail_reason[i]))
  }
}

n_pass  <- sum(ca_cell_params_all$qc_pass)
n_total <- nrow(ca_cell_params_all)
n_good  <- sum(ca_cell_params_all$qc_pass & ca_cell_params_all$signal_class == "good_paced")
n_irreg <- sum(ca_cell_params_all$qc_pass & ca_cell_params_all$signal_class == "irregular")
cat(sprintf("  -> Calcium: %d / %d cells passing QC (%.1f%%)\n",
            n_pass, n_total, 100 * n_pass / n_total))
cat(sprintf("     good_paced: %d | irregular: %d\n\n", n_good, n_irreg))


# ============================================================
# Supplementary Table 1 — well & cell counts (per group + per well)
# Always runs (no toggle). Both well_qc_gate tables are guaranteed to
# exist by this point — they are computed unconditionally above,
# regardless of run_contraction / run_calcium / run_figures.
# ============================================================

cat("Writing Supplementary Table 1...\n")
cat("────────────────────────────────────────\n")

build_supp1_group <- function(gate, exclude_wells, modality_label) {
  gate %>%
    mutate(
      is_manual    = well_uid %in% exclude_wells,
      is_failed_qc = !well_qc_pass & !is_manual
    ) %>%
    group_by(Cell_Line, Condition) %>%
    summarise(
      Modality                = modality_label,
      n_wells_imaged          = n(),
      n_wells_excluded_manual = sum(is_manual),
      n_wells_failed_qc       = sum(is_failed_qc),
      n_wells_passed          = sum(well_qc_pass),
      n_cells_segmented       = sum(n_total[well_qc_pass]),
      n_responsive_cells      = sum(n_good_paced[well_qc_pass]),
      pct_responsive          = round(
        100 * n_responsive_cells / pmax(n_cells_segmented, 1), 1),
      responsive_per_well_median = if (any(well_qc_pass))
        as.integer(median(n_good_paced[well_qc_pass])) else NA_integer_,
      responsive_per_well_min    = if (any(well_qc_pass))
        min(n_good_paced[well_qc_pass]) else NA_integer_,
      responsive_per_well_max    = if (any(well_qc_pass))
        max(n_good_paced[well_qc_pass]) else NA_integer_,
      .groups = "drop"
    ) %>%
    select(Modality, Cell_Line, Condition, everything())
}

supp_table1_groups <- bind_rows(
  build_supp1_group(ca_well_qc_gate,    ca_exclude_wells,    "Calcium"),
  build_supp1_group(contr_well_qc_gate, contr_exclude_wells, "Contraction")
) %>%
  arrange(Modality, Cell_Line, Condition)

supp_table1_wells <- bind_rows(
  ca_well_qc_gate %>%
    mutate(qc_status = case_when(
      !well_qc_pass                   ~ "failed well QC",
      well_uid %in% ca_exclude_wells  ~ "manually excluded",
      TRUE                            ~ "passed"
    )),
  contr_well_qc_gate %>%
    mutate(qc_status = case_when(
      !well_qc_pass                      ~ "failed well QC",
      well_uid %in% contr_exclude_wells  ~ "manually excluded",
      TRUE                               ~ "passed"
    ))
) %>%
  arrange(modality, Cell_Line, Condition, well_uid) %>%
  # pct_responder is calcium-only (the contraction pipeline emits no such flag),
  # so it is NA for contraction rows. Reported so the QC gate is auditable.
  select(modality, well_uid, Plate, Well, Cell_Line, Condition,
         n_total, n_good_paced, n_irregular,
         pct_good_paced, pct_irregular, pct_responder,
         qc_status, qc_fail_reason)

# Write stacked sheet to a dedicated workbook
supp1_dir <- file.path(main_dir, "Tables")
dir.create(supp1_dir, showWarnings = FALSE, recursive = TRUE)
supp1_path <- file.path(supp1_dir,
  paste0("supplementary_table_1_", analysis_date, ".xlsx"))

supp1_wb <- createWorkbook()
addWorksheet(supp1_wb, "SuppTable1")
groups_df <- relabel_cell_line_table(supp_table1_groups)
wells_df  <- relabel_cell_line_table(supp_table1_wells)

header_style <- createStyle(textDecoration = "bold",
                            fgFill = "#D9E1F2", border = "Bottom")
label_style  <- createStyle(textDecoration = "bold", fontSize = 12,
                            fgFill = "#305496", fontColour = "#FFFFFF")

writeData(supp1_wb, "SuppTable1",
          "Sub-table A — Per-group summary",
          startRow = 1, startCol = 1)
addStyle(supp1_wb, "SuppTable1", label_style, rows = 1, cols = 1)
writeData(supp1_wb, "SuppTable1", groups_df, startRow = 2, startCol = 1)
addStyle(supp1_wb, "SuppTable1", header_style,
         rows = 2, cols = seq_len(ncol(groups_df)), gridExpand = TRUE)

b_label_row <- 2 + nrow(groups_df) + 2
b_data_row  <- b_label_row + 1
writeData(supp1_wb, "SuppTable1",
          "Sub-table B — Per-well detail (all imaged wells)",
          startRow = b_label_row, startCol = 1)
addStyle(supp1_wb, "SuppTable1", label_style,
         rows = b_label_row, cols = 1)
writeData(supp1_wb, "SuppTable1", wells_df,
          startRow = b_data_row, startCol = 1)
addStyle(supp1_wb, "SuppTable1", header_style,
         rows = b_data_row, cols = seq_len(ncol(wells_df)),
         gridExpand = TRUE)

setColWidths(supp1_wb, "SuppTable1",
             cols = seq_len(max(ncol(groups_df), ncol(wells_df))),
             widths = "auto")
freezePane(supp1_wb, "SuppTable1", firstActiveRow = 3)

saveWorkbook(supp1_wb, supp1_path, overwrite = TRUE)
cat(sprintf("  Saved: %s\n\n", basename(supp1_path)))


# ============================================================
# STEPS 3 & 4: QC plots + trace galleries (optional)
# Needs the well_qc_gate tables computed above.
# ============================================================

if (run_qc_plots) {
  cat("STEP 3: Generating contraction QC plots...\n")
  cat("────────────────────────────────────────\n")
  source("3_qc_plots_contraction.R")

  cat("STEP 4: Generating calcium QC plots...\n")
  cat("────────────────────────────────────────\n")
  source("4_qc_plots_calcium.R")
} else {
  cat("STEPS 3-4: Skipped (run_qc_plots = FALSE)\n\n")
}


# ============================================================
# STEPS 5-7: Manuscript figures, tables and statistics
# ============================================================

if (run_figures) {
  cat("STEP 5: Statistics, Figures 2, 4, 5, 6 and supplementary figures...\n")
  cat("────────────────────────────────────────\n")
  source("5_plot_and_stats.R")
  cat("\nSTEP 6: Figure 3 (variance decomposition + well-average vs single-cell)...\n")
  cat("────────────────────────────────────────\n")
  source("6_figure3_variance_decomposition.R")
  cat("\nSTEP 7: Figure 1 trace panels...\n")
  cat("────────────────────────────────────────\n")
  source("7_figure1_trace_illustration.R")
} else {
  cat("STEPS 5-7: Skipped (run_figures = FALSE)\n\n")
}


# ============================================================
# Paper statistics summary
# ============================================================

cat("========================================\n")
cat("  Paper statistics summary\n")
cat("========================================\n\n")

summarise_modality <- function(params, modality_label) {
  per_well <- params %>%
    group_by(well_uid, Condition, Cell_Line) %>%
    summarise(
      n_cells      = n(),
      n_pass       = sum(qc_pass),
      n_good_paced = sum(signal_class == "good_paced"),
      n_irregular  = sum(signal_class == "irregular"),
      pct_pass     = 100 * n_pass / n_cells,
      .groups = "drop"
    )

  overall <- tibble(
    modality        = modality_label,
    n_wells         = nrow(per_well),
    n_cells_total   = sum(per_well$n_cells),
    mean_cells_well = round(mean(per_well$n_cells), 1),
    sd_cells_well   = round(sd(per_well$n_cells), 1),
    n_qc_pass       = sum(per_well$n_pass),
    pct_qc_pass     = round(100 * sum(per_well$n_pass) /
                              sum(per_well$n_cells), 1),
    mean_pct_pass_well = round(mean(per_well$pct_pass), 1),
    sd_pct_pass_well   = round(sd(per_well$pct_pass), 1),
    n_good_paced    = sum(per_well$n_good_paced),
    pct_good_paced  = round(100 * sum(per_well$n_good_paced) /
                              sum(per_well$n_cells), 1),
    n_irregular     = sum(per_well$n_irregular),
    pct_irregular   = round(100 * sum(per_well$n_irregular) /
                              sum(per_well$n_cells), 1)
  )

  by_group <- params %>%
    group_by(Condition, Cell_Line) %>%
    summarise(
      n_cells      = n(),
      n_wells      = n_distinct(well_uid),
      n_good_paced = sum(signal_class == "good_paced"),
      n_irregular  = sum(signal_class == "irregular"),
      pct_pass     = round(100 * sum(qc_pass) / n(), 1),
      .groups = "drop"
    ) %>%
    mutate(modality = modality_label)

  list(overall = overall, per_well = per_well, by_group = by_group)
}

stats_ca    <- summarise_modality(ca_cell_params_all,    "Calcium")
stats_contr <- summarise_modality(contr_cell_params_all, "Contraction")

# Print overall
overall_tbl <- bind_rows(stats_ca$overall, stats_contr$overall)
cat("--- Overall ---\n")
print(t(overall_tbl), quote = FALSE)

# Print per condition × cell line
cat("\n--- Per condition × cell line ---\n")
by_group_tbl <- bind_rows(stats_ca$by_group, stats_contr$by_group) %>%
  arrange(modality, Cell_Line, Condition)
print(as.data.frame(by_group_tbl), row.names = FALSE)

# Save
table_dir <- file.path(main_dir, "Tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(overall_tbl,
          file.path(table_dir,
                    paste0("paper_stats_overall_", analysis_date, ".csv")),
          row.names = FALSE)
write.csv(relabel_cell_line_table(by_group_tbl),
          file.path(table_dir,
                    paste0("paper_stats_by_group_", analysis_date, ".csv")),
          row.names = FALSE)
cat("\n  Saved: paper_stats_overall + paper_stats_by_group\n\n")

# ============================================================
# Done
# ============================================================

cat("========================================\n")
cat("  Pipeline complete!\n")
cat("========================================\n")
cat("\nOutputs:\n")
cat("  Figures (main):          ", pub_dir,                       "\n")
cat("  Figures (supplementary): ", expl_dir,                      "\n")
cat("  Tables:                  ", file.path(main_dir, "Tables"), "\n")
cat("  QC reports:              ", file.path(main_dir, "QC"),     "\n")
cat("\n")

# Record R and package versions used for this run
writeLines(capture.output(sessionInfo()),
           file.path(main_dir, "Tables",
                     paste0("sessionInfo_", analysis_date, ".txt")))
