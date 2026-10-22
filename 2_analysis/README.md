# 2_analysis

All pipeline scripts. Run them through `run_pipeline.R`.

| # | Script | What it does |
|---|---|---|
| — | `run_pipeline.R` | Master orchestrator. Loads packages and shared utilities, sets the four `run_*` flags, applies the manual and automatic well-level QC, writes Supplementary Table 1 and the overall cell counts, and sources steps 1–7. |
| 1 | `1_process_contraction.R` | Reads CyteSeer contraction traces + TMRM-channel BODIPY. Parses traces, detects beats, classifies cells (Otsu AC1 + beat-count cap), writes `contraction_all_cells.rds` and `contraction_results_all.rds`. |
| 2 | `2_process_calcium.R` | Reads CyteSeer calcium traces + cell/nucleus area. Extracts amplitude / CTD50 / CTD90 / rise / decay, classifies cells, writes `calcium_all_cells.rds` and `calcium_results_all.rds`. |
| 3 | `3_qc_plots_contraction.R` | Optional QC diagnostics for the contraction session (AC1 distributions, trace galleries per well). |
| 4 | `4_qc_plots_calcium.R` | Same as step 3 for the calcium session. |
| 5 | `5_plot_and_stats.R` | Mixed-effects models (`lme4` + `lmerTest`) with Tukey post-hoc tests (`emmeans`) and residual-skewness checks; Figure 2A–B; single-cell Spearman correlations (Supp. Fig. 3) and well-level calcium × contraction correlations (Supp. Fig. 4); irregular-cell GLMMs (Supp. Fig. 1); well-average two-way ANOVA vs mixed model; k-means clustering (Figs. 4–5, Supp. Table 4); within-cluster condition tests (Fig. 6); lipid-loaded vs neighbour comparison in oleate wells; `statistical_results_<date>.xlsx`. |
| 6 | `6_figure3_variance_decomposition.R` | Variance decomposition (cell line, condition, well, cell) and well-average vs single-cell panels (Figure 3). Appends the `Variance_Decomp` sheet. |
| 7 | `7_figure1_trace_illustration.R` | Eight single-cell traces per modality from one well (Figure 1 trace panels). |
| 8 | `8_figure1_parameter_schematic.R` | Parameter schematic: one typical control Ca2+ transient and contraction beat with amplitude, rise, decay, CTD50 and CTD90 marked. |

## Shared utilities (in `../0_setup/`)

| File | Contents |
|---|---|
| `colors.R` | Condition palette, cluster palettes per session, `contrast_ink()`. |
| `theme.R` | `theme_project()` (general) and `journal_theme` (figure panels). |
| `helpers.R` | Cell-line display labels, `load_plate_maps()`, `load_and_merge_cyteseer()`. |
| `contraction_trace_utils.R` | Trace parsing, beat detection, contraction + calcium parameter extraction, AC1 and Otsu threshold. |
| `cluster_helpers.R` | `run_subpop_analysis()`, `relabel_kmeans()` (fixed cluster numbering) and `plot_cluster_traces()`. |

## Working directory

Scripts use relative paths (`../1_data/`, `../3_output/`). Set the working directory to `2_analysis/` before sourcing `run_pipeline.R`.
