#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "==> Checking required files"
required=(
  LICENSE
  README.md
  CITATION.cff
  bulk/deseq2_analysis.R
  bulk/deseq2_results_NPC.txt
  bulk/deseq2_results_mESC.txt
  scRNAseq/1_QualityControl.ipynb
  scRNAseq/2_clustering.R
  scRNAseq/good_cells.csv
  4fSAMMYseq/Barplot_DiffSolGenes.R
  4fSAMMYseq/TPM_DSR.R
  reproducibility/software_versions.tsv
  reproducibility/data_accessions.tsv
  reproducibility/HOW_TO_RUN.md
  reproducibility/environments/nextflow_pins.env
  reproducibility/environments/python_xci-venv/requirements.txt
  demo/run_demo.py
)
for f in "${required[@]}"; do
  if [[ ! -e "$f" ]]; then
    echo "MISSING: $f" >&2
    exit 1
  fi
done
echo "OK"

echo "==> Demo"
python3 -m py_compile demo/run_demo.py
python3 demo/run_demo.py

echo "All checks passed."
