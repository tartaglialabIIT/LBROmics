# How to run

Scripts may contain historical absolute paths; point them at your local GEO downloads before running.

## Demo

```bash
python3 demo/run_demo.py
```

## Bulk RNA-seq (clone A8)

1. Place count matrices from GSE318892 as `bulk/raw_counts_NPC.txt` and `bulk/raw_counts_mESC.txt`, and GRCm38.98 GTF under `bulk/ref/`.
2. Optional preprocessing: `bulk/preprocessing/` (see `bulk_RNA_seq_pipeline.txt`).
3. `Rscript bulk/deseq2_analysis.R`
4. `Rscript bulk/webgestalt.R`
5. `Rscript bulk/karyoplot_mouse.R`
6. Chromosomal enrichment: `bulk/chromosomal_enrichment/Chromosome_Enrichment.ipynb`
7. Liver analyses: `bulk/liver_young_kumar/`
8. Clone B3: `bulk/clone_B3/`

## Single-cell RNA-seq

1. Cell Ranger: `scRNAseq/cellranger_count_script.sh`
2. QC: `scRNAseq/1_QualityControl.ipynb`
3. Clustering → markers → trajectories → tradeSeq → gProfiler: `2_` … `7_`
4. Cells passing QC: `scRNAseq/good_cells.csv`
5. Velocity / CellRank: `scRNAseq/RNA_velocity/`
6. In silico bulk: `scRNAseq/in_silico_bulk_DE/`, `scRNAseq/Integrated_PCA/`

## 4f-SAMMY-seq

1. Run nf-core/sammyseq using `reproducibility/environments/nextflow/cmd_to_run.sh` and pins in `nextflow_pins.env`.
2. `4fSAMMYseq/Barplot_DiffSolGenes.R`
3. `4fSAMMYseq/TPM_DSR.R`
4. `4fSAMMYseq/sammy_vs_chip.R`, `4fSAMMYseq/plot_tracks_Gviz.R`
5. Peak vs random: `4fSAMMYseq/ChIP_peaks_vs_random/`
