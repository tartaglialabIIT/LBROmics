# Live R `sessionInfo()` dumps (current machine — secondary)

**Not for manuscript version claims.** Use [`../r_sessionInfo_manuscript/`](../r_sessionInfo_manuscript/) for Methods-aligned versions.

These files were produced on 2026-09-21 by loading `library()` imports on the author laptop
(R **4.4.2** / Bioconductor **3.20** / macOS Sequoia arm64). Many package versions **differ** from the preprint Methods
(e.g. DESeq2 1.46.0 vs **1.30.1**; Seurat 5.4.0 vs **4.1.0**). See [`../manuscript_vs_current_machine.tsv`](../manuscript_vs_current_machine.tsv).

Kept only as a smoke-test / modern install aid. `destiny` failed to load (missing `smoother`).

| File | Scripts covered |
|------|-----------------|
| `bulk_RNAseq_DESeq2_sessionInfo.txt` | `bulk/deseq2_analysis.R`, clone B3 DESeq2 |
| `bulk_RNAseq_WebGestalt_sessionInfo.txt` | `bulk/webgestalt.R` |
| `bulk_RNAseq_karyoploteR_sessionInfo.txt` | karyoplot scripts |
| `bulk_RNAseq_liver_sessionInfo.txt` | liver Young/Kumar DE |
| `sammyseq_visualization_sessionInfo.txt` | `4fSAMMYseq/` viz imports |
| `scRNAseq_Seurat_stack_sessionInfo.txt` | `scRNAseq/2_`…`7_` |
| `scRNAseq_insilico_bulk_sessionInfo.txt` | in silico bulk DE |
