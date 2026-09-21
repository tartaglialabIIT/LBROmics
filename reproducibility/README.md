# Supporting materials

| File | Description |
|------|-------------|
| `data_accessions.tsv` | GEO accessions linked to repository directories |
| `software_versions.tsv` | Software versions used in the study |
| `HOW_TO_RUN.md` | Suggested execution order |
| `environments/` | Python requirements, Nextflow pins, R package inventories |

Python environment recreation (approximate):

```bash
python3.7 -m venv xci-venv
source xci-venv/bin/activate
pip install -r environments/python_xci-venv/requirements.txt
```


## Live R `sessionInfo()` dumps

Directory [`environments/r_sessionInfo_live/`](environments/r_sessionInfo_live/) contains **live** `sessionInfo()` captures obtained by loading the `library()` imports used in each workflow's R scripts on the author laptop (R 4.4.2 / Bioconductor 3.20, macOS arm64, 2026-09-21).

These are useful for reviewers installing a **modern working stack**. They are **not** bit-for-bit recreations of the manuscript environments (Methods cite older versions such as DESeq2 1.30.1 and Seurat 4.1.0). Prefer replacing with dumps from the original analysis machines when available.

Note: `destiny` failed to load here (missing CRAN dependency `smoother` for this R version); DiffusionMap-related scripts need that package resolved separately.
