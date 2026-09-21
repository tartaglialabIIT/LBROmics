#!/usr/bin/env python3
"""LBROmics demo: sanity-check deposited bulk DESeq2 result tables.

Does not re-run statistical analyses. Thresholds match the manuscript /
deseq2_analysis.R: abs(log2FoldChange) > 1 and padj < 0.01.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def read_deseq(path: Path):
    rows = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def parse_float(x: str):
    if x is None or x == "" or x.upper() == "NA":
        return float("nan")
    return float(x)


def count_degs(rows, lfc_thr: float = 1.0, padj_thr: float = 0.01):
    up = down = tested = 0
    for row in rows:
        lfc = parse_float(row.get("log2FoldChange", "NA"))
        padj = parse_float(row.get("padj", "NA"))
        if math.isnan(lfc) or math.isnan(padj):
            continue
        tested += 1
        if abs(lfc) > lfc_thr and padj < padj_thr:
            if lfc > 0:
                up += 1
            else:
                down += 1
    return tested, up, down


def write_volcano(rows, out_png: Path) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False

    xs, ys, colors = [], [], []
    for row in rows:
        lfc = parse_float(row.get("log2FoldChange", "NA"))
        padj = parse_float(row.get("padj", "NA"))
        if math.isnan(lfc) or math.isnan(padj) or padj <= 0:
            continue
        y = -math.log10(padj)
        xs.append(lfc)
        ys.append(y)
        sig = abs(lfc) > 1 and padj < 0.01
        colors.append("#F8766D" if sig and lfc > 0 else "#00BFC4" if sig else "#B0B0B0")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(xs, ys, c=colors, s=6, alpha=0.5, linewidths=0)
    ax.axvline(-1, color="grey", ls="--", lw=0.8)
    ax.axvline(1, color="grey", ls="--", lw=0.8)
    ax.axhline(-math.log10(0.01), color="grey", ls="--", lw=0.8)
    ax.set_xlabel("log2 fold change (NT-KO vs WT)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("LBROmics demo: NPC DESeq2 results")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return True


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: demo/output)",
    )
    args = parser.parse_args(argv)

    root = repo_root()
    outdir = args.outdir or (root / "demo" / "output")
    outdir.mkdir(parents=True, exist_ok=True)

    npc_path = root / "bulk" / "deseq2_results_NPC.txt"
    mesc_path = root / "bulk" / "deseq2_results_mESC.txt"
    for p in (npc_path, mesc_path):
        if not p.is_file():
            print(f"ERROR: missing required file: {p}", file=sys.stderr)
            return 1

    t0 = time.time()
    npc = read_deseq(npc_path)
    mesc = read_deseq(mesc_path)
    npc_tested, npc_up, npc_down = count_degs(npc)
    mesc_tested, mesc_up, mesc_down = count_degs(mesc)

    # Manuscript Methods (clone A8): mESC 300 up / 406 down; NPC 681 up / 828 down
    expected = {
        "mESC_up": 300,
        "mESC_down": 406,
        "NPC_up": 681,
        "NPC_down": 828,
    }

    lines = []
    lines.append("LBROmics demo summary")
    lines.append(f"NPC genes with finite LFC+padj: {npc_tested}")
    lines.append(f"NPC DEGs up (LFC>1, padj<0.01): {npc_up} (manuscript: {expected['NPC_up']})")
    lines.append(f"NPC DEGs down (LFC<-1, padj<0.01): {npc_down} (manuscript: {expected['NPC_down']})")
    lines.append(f"mESC DEGs up: {mesc_up} (manuscript: {expected['mESC_up']})")
    lines.append(f"mESC DEGs down: {mesc_down} (manuscript: {expected['mESC_down']})")

    ok = (
        npc_up == expected["NPC_up"]
        and npc_down == expected["NPC_down"]
        and mesc_up == expected["mESC_up"]
        and mesc_down == expected["mESC_down"]
    )
    lines.append("DEG count check: " + ("PASS" if ok else "FAIL"))

    volcano = outdir / "volcano_NPC.png"
    if write_volcano(npc, volcano):
        lines.append(f"Wrote {volcano.relative_to(root)}")
    else:
        lines.append("matplotlib not available; skipped volcano_NPC.png")

    elapsed = time.time() - t0
    lines.append(f"Runtime_seconds: {elapsed:.3f}")

    summary = outdir / "demo_summary.txt"
    summary.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    # Also refresh expected_output summary counts (without runtime line) for docs
    expected_path = root / "demo" / "expected_output" / "demo_summary.txt"
    expected_lines = [ln for ln in lines if not ln.startswith("Runtime_seconds")]
    expected_path.write_text("\n".join(expected_lines) + "\n")

    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
