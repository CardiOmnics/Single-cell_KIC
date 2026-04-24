# 1_data

Small input files required by the pipeline:

| File | Source | Used by |
|---|---|---|
| `Plate_Map_1.6.txt` | Manual | `1_process_contraction.R`, `2_process_calcium.R` (via `load_plate_maps()` in `0_setup/helpers.R`) |
| `Plate_Map_1.6.24.txt` | Manual | Same as above |
| `troponin.csv` | Troponin-I ELISA, one pooled sample per condition × cell line | `5_plot_and_stats.R` — Figure troponin + Bodipy |

## Raw CyteSeer data (NOT in this repository)

The raw single-cell trace exports from CyteSeer (Vala Sciences) are archived on Zenodo:

> **DOI: 10.5281/zenodo.XXXXXXX** *(placeholder — update once deposited)*

After downloading, place the Zenodo archives so the paths referenced at the top of `1_process_contraction.R` and `2_process_calcium.R` resolve:

```
1_data/
├── Calcium/
│   ├── CyteSeerResults_1.6_B02_B04/
│   ├── CyteSeerResults_1.6_B05_C04/
│   ├── CyteSeerResults_1.6_C05_E02/
│   ├── CyteSeerResults_1.6_E03_E07/
│   ├── CyteSeerResults_1.6.24_B02_B11/
│   ├── CyteSeerResults_1.6.24_C02_C11/
│   ├── CyteSeerResults_1.6.24_D02_D06/
│   ├── CyteSeerResults_1.6.24_D07_D11/
│   └── CyteSeerResults_1.6.24_E02_E11/
└── Contraction/
    ├── CyteSeerResults_2025-11-17_15-38-19_1.6/
    ├── CyteSeerResults_2025-11-17_16-56-10_1.6.24_B/
    └── CyteSeerResults_2025-12-04_12-50-13_1.6.24_CDE/
```

Each directory contains one `Primary_DataTable.csv` per imaged well.

If you only want to regenerate the publication figures, you do **not** need the raw data — the preprocessed `.rds` files in `3_output/processed/` are sufficient. Set `run_contraction = FALSE` and `run_calcium = FALSE` in `run_pipeline.R`.
