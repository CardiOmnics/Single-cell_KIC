# Single-cell functional framework for cellular state distribution profiling in hiPSC-derived cardiomyocytes

R analysis code accompanying the manuscript

> Spanjersberg T.C.F., van der Wilt C.N., *et al.* (2026).
> *Single-cell functional framework for cellular state distribution profiling in hiPSC-derived cardiomyocytes.*

Single-cell calcium, contraction, and lipid-accumulation (BODIPY) traces from hiPSC-derived cardiomyocytes are extracted from CyteSeer output, quality-controlled, classified as regular vs. irregular, and analysed with mixed-effects models, variance decomposition, and k-means subpopulation clustering.

## Repository layout

```
.
├── 0_setup/                Shared utilities (colors, theme, loader and trace helpers, clustering functions)
├── 1_data/                 Plate maps and troponin measurements (raw CyteSeer traces are archived on Zenodo — see 1_data/README.md)
├── 2_analysis/             Pipeline scripts
├── 3_output/               Written by the pipeline (processed RDS files are included; figures/tables are regenerated on run)
└── README.md               This file
```

## How to run

1. **R requirements** — R ≥ 4.4.2. The pipeline uses `tidyverse`, `ggpubr`, `broom`, `lme4`, `lmerTest`, `broom.mixed`, `emmeans`, `openxlsx`, `pracma`, `ggcorrplot`, `ggrain`, `scales`, `cluster`, `car`, `patchwork`, `cowplot`, `colorspace`. Missing packages are installed on first run.
2. **Open** `2_analysis/run_pipeline.R` with `2_analysis/` as the working directory.
3. **Configure** the four `run_*` flags near the top of `run_pipeline.R`:
   - `run_contraction` / `run_calcium` — re-process raw CyteSeer data (requires the archived raw dataset, see `1_data/README.md`). Leave `FALSE` to skip and re-use `3_output/processed/*.rds`.
   - `run_qc_plots` — regenerate QC diagnostic plots.
   - `run_figures` — regenerate all publication figures, tables, and statistical results.
4. **Source** the whole file. Outputs land under `3_output/Figures/`, `3_output/Tables/`, and `3_output/QC/`.

Re-running `run_figures = TRUE` (with the other three flags `FALSE`) only requires the six `3_output/processed/*.rds` files that ship with this repository — no raw data download needed.

## Pipeline summary

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `1_process_contraction.R` | Parse contraction traces, extract per-cell parameters + BODIPY intensity + nucleus area. |
| 2 | `2_process_calcium.R` | Parse calcium traces, extract per-cell parameters + cell/nucleus area. |
| 3 | `3_qc_plots_contraction.R` | QC diagnostic plots for the contraction session. |
| 4 | `4_qc_plots_calcium.R` | QC diagnostic plots for the calcium session. |
| 5 | `5_plot_and_stats.R` | Merge, mixed-effects models, variance decomposition (Figure 2), rainclouds, correlations, PCA, subpopulation clustering, all publication figures and tables. |

All scripts are orchestrated by `run_pipeline.R`.

## Data availability

- Plate maps and troponin-I measurements are shipped in `1_data/`.
- Processed per-cell data (six `.rds` files) ships in `3_output/processed/` so step 5 can run without the raw data.
- Raw CyteSeer single-cell trace exports (tens of GB) are archived on Zenodo — see `1_data/README.md` for the DOI.

## Citation

If you use this code please cite the paper.

## Contact

Magdalena Harakalova — m.harakalova@umcutrecht.nl
