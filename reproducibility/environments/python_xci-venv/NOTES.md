# Python environment: `xci-venv`

## Recovered locally

| Item | Value | Provenance |
|------|-------|------------|
| Python | 3.7.7 | `pyvenv.cfg` in `Cerase_single_cell/xci-venv` |
| Kernel name (QC notebook) | `xci-venv` | `scRNAseq/1_QualityControl.ipynb` |
| Canonical pin list | `requirements.txt` | `ANALYSIS/Seurat_analysis/RNA_velo_Cellrank/JAN2025/requirements.txt` |
| Full package inventory | `xci-venv.pipfreeze_from_METADATA.txt` | METADATA from local `site-packages` (88 pkgs) |

## Important caveats

1. The venv binaries still point at `/Users/jonathan/Desktop/.../xci-venv/bin/python` — `pip freeze` cannot be run live here. The METADATA freeze is the best local recovery.
2. Minor version drift vs preprint: local freeze has `scipy==1.7.3` and `statsmodels==0.13.2`; preprint Methods cite `scipy 1.7.2` and `statsmodels 0.13.1`.
3. Chromosomal enrichment notebook (`Chromosome_Enrichment.ipynb`) used kernel **`INTERACTomics-venv`** (also Python 3.7.7) and printed `statsmodels==0.13.1` — a **second** Python env, not exported here.
4. Velocity preprocessing (`kb_script.sh`) activated **`velo-venv`** on `/mnt/large/...` — not present in this OneDrive folder.

## Suggested install (approximate)

```bash
python3.7 -m venv xci-venv
source xci-venv/bin/activate
pip install -r requirements.txt
```

Prefer swapping in an author `pip freeze` from the original machine when available.
