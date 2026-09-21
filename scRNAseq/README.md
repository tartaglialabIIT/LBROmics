# scRNA-seq

| Script | Role |
|--------|------|
| `cellranger_count_script.sh` | Cell Ranger count |
| `1_QualityControl.ipynb` | Scanpy QC |
| `2_clustering.R` … `7_gProfiler_lineage2.R` | Seurat clustering through tradeSeq / gProfiler |
| `myfunctions.R` | Helper functions |
| `good_cells.csv` | QC-pass cells for clustering |
| `RNA_velocity/` | kb_python + scvelo/CellRank |
| `in_silico_bulk_DE/` | Pseudobulk DE |
| `Integrated_PCA/` | Integration with bulk RNA-seq |

Set `LBROMICS_SCRNA_ROOT` / `LBROMICS_GOOD_CELLS` as needed. Python packages: `../reproducibility/environments/python_xci-venv/`.
