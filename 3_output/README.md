# 3_output

Pipeline outputs. Apart from this README, nothing in this folder is tracked by git.

| Path | Source | Contents |
|---|---|---|
| `processed/` | Available upon request, or rebuilt by steps 1–2 from the raw data | Four `.rds` files with per-cell parameters and traces. Steps 5–7 run from these without the raw CyteSeer data. |
| `Figures/Publication/` | Pipeline | Main-figure PDFs (Figures 1–6). |
| `Figures/Exploratory/` | Pipeline | Supplementary figures, single panels and diagnostics. |
| `Tables/` | Pipeline | Supplementary Table 1, `statistical_results_<date>.xlsx`, per-cell CSVs, overall counts, `sessionInfo_<date>.txt`. |
| `QC/` | Pipeline | Normality QQ plots and (with `run_qc_plots = TRUE`) QC diagnostics. |

See the README in the repository root for which file belongs to which figure or table in the paper.

## `processed/` contents

| File | Written by | Description |
|---|---|---|
| `contraction_all_cells.rds` | Step 1 | Per-cell contraction parameters, BODIPY intensity, nucleus area, AC1 score and signal class (one row per cell). |
| `contraction_results_all.rds` | Step 1 | Parsed contraction traces per cell (used for cluster medoid traces and Figure 1). |
| `calcium_all_cells.rds` | Step 2 | Per-cell calcium parameters (amplitude, CTD50/90, rise, decay), cell and nucleus area, CyteSeer responder flag, AC1 score and signal class. |
| `calcium_results_all.rds` | Step 2 | Parsed calcium traces per cell. |

To regenerate figures and tables without the raw CyteSeer data, leave `run_contraction = FALSE` and `run_calcium = FALSE` in `run_pipeline.R`.
