# Demo: deposited bulk DE tables

This demo does **not** re-run DESeq2. It checks that the deposited differential-expression tables in `bulk/` are present and consistent with the manuscript DEG counts, and writes a simple volcano plot for NPC results.

## Why this demo?

Full manuscript reproduction needs GEO downloads, GTF annotation, and pinned R/Bioconductor environments that are not yet locked in this repository. Journal code-sharing guidelines typically expect a small runnable example; deposited DE tables are the largest scientifically relevant artefacts already in the repo.

## Requirements

- Python 3 (stdlib is enough for the numerical summary)
- Optional: `matplotlib` for PNG volcano output

## Run

From the repository root:

```bash
python3 demo/run_demo.py
```

## Expected outputs

| File | Description |
|------|-------------|
| `demo/output/demo_summary.txt` | DEG counts at \|log2FC\| > 1 and padj < 0.01 |
| `demo/output/volcano_NPC.png` | Volcano plot (if matplotlib is installed) |

Reference copy of an expected summary (counts only): `demo/expected_output/demo_summary.txt`.

## Runtime

Measured on a Linux workstation (2026-09-21): **~0.1 s** without matplotlib; **~1.7–2 s** when writing `volcano_NPC.png`.
