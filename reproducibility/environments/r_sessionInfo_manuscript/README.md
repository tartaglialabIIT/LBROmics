# Manuscript-aligned R environment inventories

**Canonical for reproducing the paper.** Versions are taken from Fiorentino et al. bioRxiv Methods
([doi:10.64898/2026.03.30.714681](https://doi.org/10.64898/2026.03.30.714681)), plus script comments where Methods omit a version.

These files are **not** live `sessionInfo()` dumps from an original analysis machine.
They are Methods-faithful inventories in `sessionInfo`-like form so reviewers can see which stack belongs to which workflow.

| File | Analysis stack | R (as stated / implied) |
|------|----------------|-------------------------|
| `bulk_RNAseq_sessionInfo_manuscript.txt` | DESeq2 / WebGestaltR / karyoploteR / liver DE | Not stated; DESeq2 1.30.1 ≈ Bioconductor 3.12 / R 4.0.x |
| `scRNAseq_sessionInfo_manuscript.txt` | Seurat → destiny / slingshot / condiments / tradeSeq | **R 4.2.2** (script comment) |
| `sammyseq_viz_sessionInfo_manuscript.txt` | Gviz / rtracklayer visualization | **R 4.4.2** (Methods) |

Also see:

- `packages_by_workflow.tsv` — flat table of Package / Version / Workflow / Source
- `../manuscript_vs_current_machine.tsv` — where a live laptop capture differs from Methods
- `../r_sessionInfo_live/` — optional modern-machine smoke-test dumps (**do not use for manuscript claims**)

When original `sessionInfo()` exports become available, replace these summaries with the real dumps.
