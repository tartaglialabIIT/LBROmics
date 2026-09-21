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
