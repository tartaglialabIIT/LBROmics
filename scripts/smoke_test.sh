#!/usr/bin/env bash
# Smoke tests for LBROmics reproducibility packaging.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "==> Checking required tracked files"
required=(
  LICENSE
  README.md
  CITATION.cff
  bulk/deseq2_analysis.R
  bulk/deseq2_results_NPC.txt
  bulk/deseq2_results_mESC.txt
  bulk/genes_for_karyo.csv
  scRNAseq/1_QualityControl.ipynb
  scRNAseq/2_clustering.R
  4fSAMMYseq/Barplot_DiffSolGenes.R
  4fSAMMYseq/TPM_DSR.R
  reproducibility/software_versions.tsv
  reproducibility/data_accessions.tsv
  reproducibility/figure_to_code.tsv
  demo/run_demo.py
)
for f in "${required[@]}"; do
  if [[ ! -e "$f" ]]; then
    echo "MISSING: $f" >&2
    exit 1
  fi
done
echo "OK: required files present"

echo "==> Python syntax check (demo)"
python3 -m py_compile demo/run_demo.py

echo "==> Running demo"
python3 demo/run_demo.py

echo "==> Checking for leftover absolute author paths in key scripts"
# Fail if old absolute paths remain in portability-patched entrypoints
if grep -nE '/Users/jonathan|/mnt/large/jfiorentino|/Users/jfiorentino' \
  bulk/karyoplot_mouse.R \
  bulk/webgestalt.R \
  scRNAseq/2_clustering.R \
  scRNAseq/3_cluster_markers.R \
  scRNAseq/4_diffmap.R \
  scRNAseq/5_trajectory_inference_and_condiments.R \
  scRNAseq/6_run_tradeseq.R \
  scRNAseq/7_gProfiler_lineage2.R \
  4fSAMMYseq/Barplot_DiffSolGenes.R \
  4fSAMMYseq/TPM_DSR.R; then
  echo "ERROR: absolute author paths still present in patched scripts" >&2
  exit 1
fi
echo "OK: patched scripts have no hard-coded author absolute paths"

echo "All smoke tests passed."
