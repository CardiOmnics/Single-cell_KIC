# 3_output

Pipeline outputs. Most subdirectories are **regenerated** on every run; only `processed/` is shipped with the repository.

| Path | Ships in repo? | Contents |
|---|---|---|
| `processed/` | **Yes** | Six `.rds` files holding per-cell parameters and traces. Allow step 5 to run without the raw CyteSeer data. |
| `Figures/` | No | Publication PDFs written by step 5. |
| `Tables/` | No | CSV / XLSX statistics tables written by step 5. |
| `QC/` | No | QC diagnostic plots written by steps 3–4. |

## `processed/` contents

| File | Written by | Description |
|---|---|---|
| `contraction_all_cells.rds` | Step 1 | Parsed contraction traces (one row per time point per cell). |
| `contraction_results_all.rds` | Step 1 | Per-cell contraction parameters + BODIPY intensity + nucleus area + QC flags. |
| `calcium_all_cells.rds` | Step 2 | Parsed, normalised calcium traces. |
| `calcium_results_all.rds` | Step 2 | Per-cell calcium parameters (amplitude, CTD50/90, rise, decay) + QC flags. |
| `morphology_bodipy_all_cells.rds` | Step 1 | Per-cell morphology + BODIPY measurements from the contraction session. |
| `pacing_fidelity_scores.rds` | Step 1 | Per-cell pacing-fidelity score (AC1) used in the QC classification. |

To regenerate figures/tables without touching the raw CyteSeer data, set `run_contraction = FALSE` and `run_calcium = FALSE` in `run_pipeline.R` and leave `run_figures = TRUE`.
