# How to run

Scripts often contain historical absolute paths; point them at your local copies of GEO data before running.

## Demo (no download required)

```bash
python3 demo/run_demo.py
```

## Bulk RNA-seq (clone A8)

1. Place count matrices from GSE318892 as `bulk/raw_counts_NPC.txt` and `bulk/raw_counts_mESC.txt`, and GRCm38.98 GTF under `bulk/ref/`.
2. Optional preprocessing from FASTQs: `bulk/preprocessing/` (see `bulk_RNA_seq_pipeline.txt`).
3. `Rscript bulk/deseq2_analysis.R`
4. `Rscript bulk/webgestalt.R`
5. `Rscript bulk/karyoplot_mouse.R`
6. Chromosomal enrichment: `bulk/chromosomal_enrichment/Chromosome_Enrichment.ipynb`
7. Liver analyses: `bulk/liver_young_kumar/`
8. Clone B3 companions: `bulk/clone_B3/`

## Single-cell RNA-seq

1. Cell Ranger: `scRNAseq/cellranger_count_script.sh`
2. QC: `scRNAseq/1_QualityControl.ipynb` (Python env: `reproducibility/environments/python_xci-venv/`)
3. Clustering → markers → trajectories → tradeSeq → gProfiler: `2_` … `7_`
4. Use deposited `scRNAseq/good_cells.csv` (or set `LBROMICS_GOOD_CELLS`)
5. Velocity / CellRank: `scRNAseq/RNA_velocity/`
6. In silico bulk integration: `scRNAseq/in_silico_bulk_DE/`, `scRNAseq/Integrated_PCA/`

## 4f-SAMMY-seq

1. Run nf-core/sammyseq using pins in `reproducibility/environments/nextflow_pins.env` and the samplesheet under `4fSAMMYseq/pipeline_launch/`.
2. Differential solubility gene counts: `4fSAMMYseq/Barplot_DiffSolGenes.R`
3. TPM in DSRs: `4fSAMMYseq/TPM_DSR.R`
4. Compartments: `4fSAMMYseq/compartments_analysis.R`
5. ChIP comparison / Gviz tracks: `4fSAMMYseq/sammy_vs_chip.R`, `4fSAMMYseq/plot_tracks_Gviz.R`
6. Peak vs random enrichment: `4fSAMMYseq/ChIP_peaks_vs_random/`
