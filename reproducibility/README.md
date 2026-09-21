# Reproducibility materials

| File | Purpose |
|------|---------|
| `data_accessions.tsv` | GEO / public data mapping |
| `software_versions.tsv` | Software inventory with provenance |
| `figure_to_code.tsv` | Manuscript figure ↔ script map |
| `execution_order.md` | Step-by-step run order |
| `RECOVERED_CODE.md` | Map of recovered analysis scripts |
| `environments/` | Python/Nextflow pins and R export notes |

## Environments

| Workflow | What is deposited | Notes |
|----------|-------------------|-------|
| scRNA / CellRank Python | `environments/python_xci-venv/` (`requirements.txt` + METADATA inventory); Python **3.7.7** | Prefer live `pip freeze` from original machine if available |
| 4f-SAMMY Nextflow | `environments/nextflow_pins.env` — NF **24.10.4**, sammyseq **`fa6f6ffeb3`** (`dev`), docker profile | Confirm compartments branch commit on server |
| R stacks | Preprint versions in `software_versions.tsv` / `*_r_packages.txt` | Author `sessionInfo()` still recommended (`R_sessionInfo_TODO.md`) |
| Chromosomal enrichment Python | Separate historical `INTERACTomics-venv` | See notebook kernel metadata |
| kb_python | Historical `velo-venv` on analysis server | Export still useful if available |
