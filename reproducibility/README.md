# Supporting materials

| File | Description |
|------|-------------|
| `data_accessions.tsv` | GEO accessions linked to repository directories |
| `software_versions.tsv` | Software versions used in the study (Methods) |
| `HOW_TO_RUN.md` | Suggested execution order |
| `environments/` | Python requirements, Nextflow pins, R inventories |

## R environments (manuscript-aligned)

**Use these for paper-faithful dependency claims:**

- [`environments/r_sessionInfo_manuscript/`](environments/r_sessionInfo_manuscript/) — Methods versions per workflow (bulk / scRNA-seq / SAMMY viz)
- [`environments/manuscript_vs_current_machine.tsv`](environments/manuscript_vs_current_machine.tsv) — differences vs a 2026-09-21 laptop smoke-test
- [`environments/bulk_r_packages.txt`](environments/bulk_r_packages.txt), [`environments/scrnaseq_r_packages.txt`](environments/scrnaseq_r_packages.txt)

The preprint used **multiple** R stacks (older Bioconductor for DESeq2 1.30.1; R 4.2.2 for Seurat/trajectories; R 4.4.2 for SAMMY Gviz). A single lockfile would be misleading.

Optional secondary dumps from loading imports on a current laptop (R 4.4.2 / Bioc 3.20) are in [`environments/r_sessionInfo_live/`](environments/r_sessionInfo_live/) — **do not cite those versions as manuscript Methods**.

## Python

```bash
python3.7 -m venv xci-venv
source xci-venv/bin/activate
pip install -r environments/python_xci-venv/requirements.txt
```

Scanpy/scvelo/cellrank pins also follow manuscript Methods (see `software_versions.tsv`).
