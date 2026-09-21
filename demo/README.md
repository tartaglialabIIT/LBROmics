# Demo: deposited bulk DE tables

This demo checks that the deposited differential-expression tables in `bulk/` are present and consistent with the manuscript DEG counts, and optionally writes a volcano plot for NPC results.

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

Reference summary: `demo/expected_output/demo_summary.txt`.

## Runtime

Typically ~0.1 s without matplotlib, ~2 s with volcano PNG.
