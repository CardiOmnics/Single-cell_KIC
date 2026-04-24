# ============================================================
# Master pipeline: Single-cell calcium, contraction & Bodipy
# Run from top to bottom — no interactive prompts.
# To exclude specific wells, edit the vectors below before running.
# ============================================================

rm(list = ls())
set.seed(123)

# Use cloud CRAN mirror so install.packages() never prompts for selection
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(broom)
  library(lme4)
  library(emmeans)
  library(openxlsx)
  library(lmerTest)
  library(broom.mixed)
  library(pracma)
})

# Run this script with 2_analysis/ as the working directory
# (open the .Rproj at the repo root, or setwd() to 2_analysis/ before sourcing).
# All paths below are relative to 2_analysis/.

# Global settings
main_dir <- "../3_output"
analysis_date <- format(Sys.Date(), "%Y-%m-%d")

# Load project-wide settings
source("../0_setup/colors.R")
source("../0_setup/theme.R")
source("../0_setup/helpers.R")
source("../0_setup/contraction_trace_utils.R")
source("../0_setup/cluster_helpers.R")

theme_set(theme_project())

cat("========================================\n")
cat("  Single-cell analysis pipeline\n")
cat("  Date:", analysis_date, "\n")
cat("========================================\n\n")

# ============================================================
# Steps to run — set to FALSE to skip
# ============================================================
run_contraction <- FALSE
run_calcium     <- FALSE
run_qc_plots    <- FALSE   # regenerate QC plots + trace galleries from saved RDS
run_figures     <- TRUE

# ============================================================
# Well exclusions — edit these before running if needed
# Format: "plate_well", e.g. "1.6_B02" or "1.6.24_D02"
# Currently empty: all exclusions handled by automatic well-level QC (min_good_paced).
# ============================================================
contr_exclude_wells <- c(
  "1.6.24_B04",  # beating at ~2 Hz (double pacing freq); no cells follow 1 Hz pacing
  "1.6_B03"      # beating at ~2 Hz (double pacing freq); no cells follow 1 Hz pacing
)
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
  rds <- file.path(main_dir, "processed", "contraction_all_cells.rds")
  if (!file.exists(rds)) stop("No saved contraction RDS found. Run step 1 first.")
  contr_cell_params_all <<- readRDS(rds)
  cat("  Loaded:", rds, "\n")
  results_rds <- file.path(main_dir, "processed", "contraction_results_all.rds")
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
      n_manually_excluded == n_total ~ "manually excluded",
      !well_qc_pass ~ sprintf("n_good_paced = %d (< %d)",
                               n_good_paced, min_good_paced),
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
  rds <- file.path(main_dir, "processed", "calcium_all_cells.rds")
  if (!file.exists(rds)) stop("No saved calcium RDS found. Run step 2 first.")
  ca_cell_params_all <<- readRDS(rds)
  cat("  Loaded:", rds, "\n")
  results_rds <- file.path(main_dir, "processed", "calcium_results_all.rds")
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

# Well-level QC: exclude wells with fewer than min_good_paced good_paced cells
ca_well_qc_gate <- ca_cell_params_all %>%
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
      n_manually_excluded == n_total ~ "manually excluded",
      !well_qc_pass ~ sprintf("n_good_paced = %d (< %d)",
                               n_good_paced, min_good_paced),
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
  cat(sprintf("  -> Well-level QC: excluded %d calcium wells with < %d good_paced cells, removed %d cells\n",
              nrow(ca_failed_wells), min_good_paced, n_before - n_after))
  cat("     Failed wells:", paste(ca_failed_wells$well_uid, collapse = ", "), "\n")
}

n_pass  <- sum(ca_cell_params_all$qc_pass)
n_total <- nrow(ca_cell_params_all)
n_good  <- sum(ca_cell_params_all$qc_pass & ca_cell_params_all$signal_class == "good_paced")
n_irreg <- sum(ca_cell_params_all$qc_pass & ca_cell_params_all$signal_class == "irregular")
cat(sprintf("  -> Calcium: %d / %d cells passing QC (%.1f%%)\n",
            n_pass, n_total, 100 * n_pass / n_total))
cat(sprintf("     good_paced: %d | irregular: %d\n\n", n_good, n_irreg))


# ============================================================
# STEPS 3 & 4: QC plots + trace galleries
# Always regenerated when run_qc_plots = TRUE, because the
# well QC gallery requires the well_qc_gate tables computed above.
# ============================================================

if (run_qc_plots) {
  cat("STEP 3: Generating contraction QC plots...\n")
  cat("────────────────────────────────────────\n")
  source("3_qc_plots_contraction.R")

  cat("STEP 4: Generating calcium QC plots...\n")
  cat("────────────────────────────────────────\n")
  source("4_qc_plots_calcium.R")
} else {
  cat("STEP 3/4: Skipped (run_qc_plots = FALSE)\n\n")
}


# ============================================================
# STEP 5: Generate publication figures and statistics
#         (includes Figure 2 variance decomposition + graphical-abstract
#          trace illustration as inline sections at the end of the script)
# ============================================================

if (run_figures) {
  cat("STEP 5: Generating publication figures and statistics...\n")
  cat("────────────────────────────────────────\n")
  source("5_plot_and_stats.R")
} else {
  cat("STEP 5: Skipped (run_figures = FALSE)\n\n")
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
cat("  QC reports:       3_output/QC/\n")
cat("  Figures:          3_output/Figures/\n")
cat("  Tables:           3_output/Tables/\n")
cat("  Processed data:   3_output/processed/\n")
cat("\n")
