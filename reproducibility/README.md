# Supporting materials

| File | Description |
|------|-------------|
| `data_accessions.tsv` | GEO accessions linked to repository directories |
| `software_versions.tsv` | Software versions from the manuscript Methods |
| `HOW_TO_RUN.md` | Suggested execution order |
| `environments/sessionInfo/` | Per-workflow version lists matching Methods |
| `environments/` | R/Python package pins and Nextflow launch command |

Python (approximate):

```bash
python3.7 -m venv xci-venv
source xci-venv/bin/activate
pip install -r environments/python_xci-venv/requirements.txt
```
