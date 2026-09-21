# Reproducibility materials

This directory supports Nature-style code sharing for the LBROmics analysis repository. It documents **what can be verified from the repository and preprint**, and clearly marks what still requires author-exported environments.

| File | Purpose |
|------|---------|
| `data_accessions.tsv` | GEO / public data mapping |
| `software_versions.tsv` | Software inventory with provenance |
| `figure_to_code.tsv` | Manuscript figure ↔ script map (partial; TODOs marked) |
| `environments/` | Package inventories and environment skeletons (**not** fabricated lockfiles) |
| `execution_order.md` | Step-by-step run order per workflow |

## Environment recreation status

| Workflow | Lockfile in repo? | Versions recoverable? | Action for authors |
|----------|-------------------|----------------------|--------------------|
| bulk RNA-seq (R) | No | Partial (preprint) | Export `sessionInfo()` + package versions from the DESeq2 machine |
| scRNA-seq Python QC | No | Python 3.7.7 + Scanpy 1.9.1 (notebook/preprint) | `pip freeze` from `xci-venv` |
| scRNA-seq R (Seurat → tradeSeq) | No | Partial (preprint + R 4.2.2 comment) | `sessionInfo()` from Seurat machine |
| 4f-SAMMY-seq nf-core | No | Unpinned `dev` in preprint | Pin commit + Nextflow version + profile |
| ChIP-seq nf-core | No | Not specified | Pin nf-core/chipseq release/commit |

Do **not** treat preprint version lists as a substitute for `renv.lock` / conda lockfiles when claiming bit-for-bit reproducibility.
