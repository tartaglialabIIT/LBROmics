# Execution order

Paths assume the repository root unless noted. Scripts often contain historical absolute paths — adapt them to your machine or GEO downloads.

## A. Bulk RNA-seq

### A0. Preprocessing (optional; raw FASTQs)

Scripts under `bulk/preprocessing/` (cutadapt, FastQC, rRNA depletion, STAR). See `bulk/preprocessing/bulk_RNA_seq_pipeline.txt` and shell helpers. Genome: Ensembl GRCm38.98.

### A1. Clone A8 DE + figures (`bulk/`)

1. Download counts from [GSE318892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318892) as `bulk/raw_counts_NPC.txt` and `bulk/raw_counts_mESC.txt`.
2. Provide `bulk/ref/Mus_musculus.GRCm38.98.gtf`.
3. `Rscript bulk/deseq2_analysis.R`
4. `Rscript bulk/webgestalt.R`
5. `Rscript bulk/karyoplot_mouse.R`

### A2. Chromosomal enrichment

`bulk/chromosomal_enrichment/Chromosome_Enrichment.ipynb` (Python; historically `INTERACTomics-venv`).

### A3. Liver Young / Kumar

Scripts in `bulk/liver_young_kumar/` using [GSE165447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165447) / [GSE324396](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324396).

### A4. Clone B3 companions

`bulk/clone_B3/` — review alongside A8 scripts before re-running.

---

## B. Single-cell RNA-seq (`scRNAseq/`)

1. Cell Ranger: `scRNAseq/cellranger_count_script.sh` (Cell Ranger 6.1.2 / GRCm38).
2. QC: `1_QualityControl.ipynb` with `xci-venv` (Python 3.7.7; see `reproducibility/environments/python_xci-venv/`).
3. Clustering: `2_clustering.R` — uses deposited `good_cells.csv` (or `LBROMICS_GOOD_CELLS`).
4. Markers → diffusion map → slingshot/condiments → tradeSeq → gProfiler: scripts `3_`…`7_`.
5. RNA velocity / CellRank: `RNA_velocity/` (`kb_script.sh`, `RNAVeloCellrank.ipynb`).
6. In silico bulk / integrated PCA: `in_silico_bulk_DE/`, `Integrated_PCA/`.

Set `LBROMICS_SCRNA_ROOT` / `LBROMICS_SCRNA_WT_10X` / `LBROMICS_SCRNA_MUT_10X` as needed.

---

## C. 4f-SAMMY-seq (`4fSAMMYseq/`)

### C0. nf-core/sammyseq

See `reproducibility/environments/nextflow_pins.env` and `reproducibility/environments/nextflow/cmd_to_run.sh`.  
Observed: Nextflow **24.10.4**, nf-core/sammyseq **`dev` @ `fa6f6ffeb3`**, `-profile docker`.  
Samplesheets: `4fSAMMYseq/pipeline_launch/`.  
Confirm `compartments_subworkflow` commit on the analysis server before claiming bit-for-bit identity.

### C1. Downstream visualization

| Script area | Role |
|-------------|------|
| `Barplot_DiffSolGenes.R`, `TPM_DSR.R` | Primary DiffSol / TPM–DSR plots |
| `downstream_viz/` | ChIP comparison, Gviz, additional DiffSol plots |
| `compartments_analysis.R` | Compartment plots |
| `ChIP_peaks_vs_random/` | Peak vs random enrichment |

Public ChIP tracks: [GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107). nf-core/chipseq pin still to be confirmed by authors.

---

## D. Demo (no GEO download)

```bash
python3 demo/run_demo.py
# or: make demo
```
