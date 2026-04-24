# 2_analysis

All pipeline scripts. Run them through `run_pipeline.R`.

| # | Script | What it does |
|---|---|---|
| — | `run_pipeline.R` | Master orchestrator. Sets packages, globals, and the four `run_*` flags; sources steps 1–5 in order. |
| 1 | `1_process_contraction.R` | Reads CyteSeer contraction traces + TMRM-channel BODIPY. Parses traces, detects beats, classifies cells (Otsu AC1 + beat-count cap), writes `contraction_all_cells.rds` and `contraction_results_all.rds`. |
| 2 | `2_process_calcium.R` | Reads CyteSeer calcium traces + cell/nucleus area. Parses normalised traces, extracts amplitude / CTD50 / CTD90 / rise / decay, classifies cells, writes `calcium_all_cells.rds` and `calcium_results_all.rds`. |
| 3 | `3_qc_plots_contraction.R` | QC diagnostic plots for the contraction session (per-well galleries, AC1 thresholds, beat-count caps). |
| 4 | `4_qc_plots_calcium.R` | Same as step 3 for the calcium session. |
| 5 | `5_plot_and_stats.R` | Merge sessions, mixed-effects models (`lme4` + `lmerTest`), pairwise Tukey (`emmeans`), raincloud plots, correlation panels, PCA, k-means subpopulation clustering, variance decomposition (Figure 2), well-average vs single-cell comparison (Figure 2), troponin + BODIPY combined panel, graphical-abstract trace illustration. Writes all publication PDFs and tables. |

## Shared utilities (sourced from `../0_setup/` by `run_pipeline.R`)

| File | Contents |
|---|---|
| `colors.R` | `group_colors`, `group_levels` — condition palette and factor order. |
| `theme.R` | `theme_project()` (general) and `journal_theme` (Nature-style). |
| `helpers.R` | `read_cyteseer_singlecell()`, `load_plate_maps()`, `load_and_merge_cyteseer()`. |
| `contraction_trace_utils.R` | Trace parsing, preprocessing (smoothing + detrending), beat detection, contraction + calcium parameter extraction. |
| `cluster_helpers.R` | `run_subpop_analysis()` and `plot_cluster_traces()` for the k-means clustering section of step 5. |

## Working directory

Scripts use relative paths (`../1_data/`, `../3_output/`). Open `run_pipeline.R` with `2_analysis/` as the working directory before sourcing.
