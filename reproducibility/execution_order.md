# Execution order

Paths below assume the repository root is the working directory unless noted. Replace `TODO` placeholders with GEO downloads and author machine paths.

## A. Bulk RNA-seq (clone A8) — scripts in `bulk/`

**Prerequisites**

1. Download processed/raw count matrices from GEO [GSE318892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318892).
2. Place count tables as `bulk/raw_counts_NPC.txt` and `bulk/raw_counts_mESC.txt` (first column gene_id, remaining columns samples; 12 samples expected — see header comments in `deseq2_analysis.R`).
3. Download Ensembl `Mus_musculus.GRCm38.98.gtf` into `bulk/ref/` (or set paths in the script).
4. Optional: gene position table `Mus_musculus.GRCm38.98_gen_pos.txt` for chromosome-annotated supplementary tables.

**Steps**

| Step | Command | Outputs |
|------|---------|---------|
| 1 | `cd bulk && Rscript deseq2_analysis.R` | `bulk/analysis_results/` (DE tables, PCA, heatmaps). Script `setwd`s into `analysis_results`. |
| 2 | `cd bulk && Rscript webgestalt.R` | WebGestalt ORA/GSEA under `bulk/analysis_results/webgestalt/` (after portability patch). Requires network access to WebGestalt resources. |
| 3 | `cd bulk && Rscript karyoplot_mouse.R` | `mm10_DE_Xchr_sig_genes.pdf` in `bulk/`. Uses deposited `deseq2_results_NPC.txt` + `genes_for_karyo.csv`. Needs Ensembl/biomaRt access. |

**Note:** STAR/cutadapt/Bowtie2 preprocessing and chromosomal-enrichment Python code are described in the preprint but are **not** in this repository.

**Approximate resources:** DESeq2 on 12 samples is typically feasible on a laptop (TODO: author-measured runtime/RAM). WebGestaltR uses `nThreads=20` in the script.

---

## B. Single-cell RNA-seq — `scRNAseq/`

**Prerequisites**

1. Download 10x matrices from GEO [GSE318871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318871) into directories matching `raw_data/s1_DR1` (WT) and `raw_data/s2_A8` (mutant), or set `LBROMICS_SCRNA_ROOT`.
2. Provide `good_cells.csv` (QC pass list). The tracked `good_cells.numbers` file is an Apple Numbers document and is **not** a drop-in replacement.

**Steps**

| Step | Script | Notes |
|------|--------|-------|
| 1 | `1_QualityControl.ipynb` | Set `main_folder`. Writes `good_cells.csv`. Python 3.7.7 / Scanpy 1.9.1 (preprint/notebook). |
| 2 | `2_clustering.R` | Integration, clustering, UMAPs; writes `RDS_objects/sub_data_combined_and_clustered.rds`. |
| 3 | `3_cluster_markers.R` | Markers, violins, heatmaps, GO via `myfunctions.R`. |
| 4 | `4_diffmap.R` | destiny DiffusionMap; writes RDS with `dm` reduction. |
| 5 | `5_trajectory_inference_and_condiments.R` | slingshot + condiments; comment states **R 4.2.2**. |
| 6 | `6_run_tradeseq.R` | `fitGAM` with 24 workers by default — reduce `BPPARAM$workers` on smaller machines. Expects `curves.RData`, `counts.Rdata`, `myconditions.Rdata`. |
| 7 | `7_gProfiler_lineage2.R` | Enrichment for lineage-2 conditionTest genes. |

**Missing vs Methods:** CellRank/scvelo/kb_python and limma in-silico bulk integration scripts are not deposited.

---

## C. 4f-SAMMY-seq — `4fSAMMYseq/`

**Upstream (not in repo)**

1. Run [nf-core/sammyseq](https://nf-co.re/sammyseq) on GEO [GSE318873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318873) with `--differential_solubility` and compartmentalization options as in Methods.
2. **Authors must pin** the exact git commit (Methods cite `/dev/` and `compartments_subworkflow`).
3. Reprocess Bonev et al. ChIP-seq ([GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107)) with nf-core/chipseq (version TODO).

**Downstream visualization (in repo)**

| Step | Script | Required inputs |
|------|--------|-----------------|
| 1 | `Barplot_DiffSolGenes.R` | Folder of gene-list files from differential solubility (`*S2SvsS3*`) |
| 2 | `TPM_DSR.R` | `S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds`, `mm10.gtf`, `salmon.merged.gene_tpm.tsv` |

Set `LBROMICS_SAMMY_ROOT` or pass paths after the portability patch (see script headers).

---

## D. Demo (no GEO download)

From repository root:

```bash
python3 demo/run_demo.py
# or
make demo
```

See `demo/README.md`.
