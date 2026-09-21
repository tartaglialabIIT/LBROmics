# Live R `sessionInfo()` dumps

Captured by loading `library()` imports used in repository R scripts (scripts were **not** fully executed).

| File | Scripts covered |
|------|-----------------|
| `bulk_RNAseq_DESeq2_sessionInfo.txt` | `bulk/deseq2_analysis.R`, `bulk/clone_B3/deseq2_analysis_cloneB3_context.R` |
| `bulk_RNAseq_WebGestalt_sessionInfo.txt` | `bulk/webgestalt.R` |
| `bulk_RNAseq_karyoploteR_sessionInfo.txt` | `bulk/karyoplot_mouse.R`, `bulk/clone_B3/karyoplot_mouseB3.R` |
| `bulk_RNAseq_liver_sessionInfo.txt` | `bulk/liver_young_kumar/DE_analysis_young_kumar.R` |
| `sammyseq_visualization_sessionInfo.txt` | `4fSAMMYseq/*.R` (Gviz/TPM/compartments/ChIP viz imports) |
| `scRNAseq_Seurat_stack_sessionInfo.txt` | `scRNAseq/2_`…`7_`, `myfunctions.R` |
| `scRNAseq_insilico_bulk_sessionInfo.txt` | `scRNAseq/in_silico_bulk_DE/DE_in_silico_bulk.R` |

Companion `*_packages.csv` files list Package/Version for requested packages that were installed.

## Important caveats

1. These dumps are from the **author laptop at capture time** (see headers: R version, Bioconductor, host OS).
2. They are **not** bit-for-bit recreations of the manuscript environments (Methods cite older stacks such as DESeq2 1.30.1 / Seurat 4.1.0 / R 4.2.2 for trajectories).
3. Prefer replacing with `sessionInfo()` from the **original analysis machines** when available.
4. Re-generate with:

```bash
Rscript reproducibility/environments/r_sessionInfo_live/capture_workflow_sessionInfo.R \
  WORKFLOW_NAME outfile.txt packages.csv pkg1 pkg2 ...
```
