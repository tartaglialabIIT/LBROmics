# scRNA-seq scripts

Numbered workflow: `1_QualityControl.ipynb` → `2_clustering.R` … → `7_gProfiler_lineage2.R`.

Additional recovered analyses: `RNA_velocity/`, `in_silico_bulk_DE/`, `Integrated_PCA/`, `cellranger_count_script.sh`.

## Cell QC list

- Use **`good_cells.csv`** (deposited; 3560 cells) with `2_clustering.R`.
- `good_cells.numbers` is a legacy Apple Numbers file and is not read by the R scripts.
- Override path with `LBROMICS_GOOD_CELLS` if needed.

Python environment: `../reproducibility/environments/python_xci-venv/`.
