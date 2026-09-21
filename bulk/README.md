# Bulk RNA-seq

| Path | Role |
|------|------|
| `preprocessing/` | cutadapt, FastQC, Bowtie2 rRNA filter, STAR |
| `deseq2_analysis.R` | DESeq2 DE + QC/heatmaps (clone A8) |
| `webgestalt.R` | ORA / GSEA |
| `karyoplot_mouse.R` | X-chromosome karyoplot |
| `deseq2_results_*.txt`, `genes_for_karyo.csv` | Deposited results used by demo/karyoplot |
| `chromosomal_enrichment/` | Chromosome enrichment of DEGs |
| `liver_young_kumar/` | Liver DE (Young et al. + new samples) |
| `clone_B3/` | Clone B3 DE / karyoplot companions |
