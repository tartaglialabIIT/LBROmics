# LBROmics

Analysis and visualization code accompanying:

> Fiorentino et al. *LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization.*  
> bioRxiv [doi:10.64898/2026.03.30.714681](https://doi.org/10.64898/2026.03.30.714681)

This repository is intended for **reproducibility and reviewer usability**, not as a general-purpose software package. It contains custom R and Python workflows for bulk RNA-seq, single-cell RNA-seq, and downstream 4f-SAMMY-seq visualization.

**License:** MIT  
**Repository:** https://github.com/tartaglialabIIT/LBROmics  
**Reproducibility audit:** [`REPRODUCIBILITY_AUDIT.md`](REPRODUCIBILITY_AUDIT.md)  
**Reproducibility tables:** [`reproducibility/`](reproducibility/)

---

## Overview

The nuclear lamina protein LBR is dissected genetically to separate nucleoplasmic-domain functions from sterol-reductase activity. The computational analyses here support differential expression, single-cell trajectory inference, and chromatin-solubility visualization described in the manuscript.

**Scope note:** Upstream Nextflow/nf-core runs (4f-SAMMY-seq, ChIP-seq reprocessing), CellRank/scvelo, liver DE, and chromosomal-enrichment Python code are described in the preprint Methods but are **not fully deposited** here. See the audit for gaps.

---

## Repository structure

| Path | Contents |
|------|----------|
| `bulk/` | DESeq2 differential expression, QC plots/heatmaps, WebGestaltR ORA/GSEA, X-chromosome karyoplot; deposited DE result tables |
| `scRNAseq/` | Scanpy QC notebook → Seurat clustering/markers → destiny/slingshot/condiments → tradeSeq → gProfiler |
| `4fSAMMYseq/` | Downstream visualization of differential solubility gene counts and TPM in DSRs (not the nf-core pipeline itself) |
| `demo/` | Small runnable demo using deposited DE tables (no GEO download) |
| `reproducibility/` | Accessions, software inventory, figure map, execution order, environment skeletons |
| `scripts/` | Smoke tests and helpers |
| `REPRODUCIBILITY_AUDIT.md` | Reproducibility gaps + author information still required |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | MIT |

---

## Data availability

| Analysis | Biological system | GEO accession | Raw/processed | Repo directory |
|----------|-------------------|---------------|---------------|----------------|
| Bulk RNA-seq (clone A8) | Female XX mESC / day-5 NPC | [GSE318892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318892) | raw + processed | `bulk/` |
| Bulk RNA-seq (clone B3) | Female XX mESC / day-5 NPC | [GSE318895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318895) | raw + processed | scripts TODO |
| scRNA-seq | WT vs Lbr NT-KO day-5 NPC | [GSE318871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318871) | 10x | `scRNAseq/` |
| 4f-SAMMY-seq | ESC / NPC clone B3 | [GSE318873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318873) | raw + processed | `4fSAMMYseq/` (viz only) |
| Bulk RNA-seq liver (new) | Male/female liver | [GSE324396](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324396) | raw + processed | scripts TODO |
| Bulk RNA-seq liver (prior) | Young et al. subset | [GSE165447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165447) | public | scripts TODO |
| ChIP-seq (public) | Bonev et al. NPCs | [GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107) | BigWig/peaks | scripts TODO |

Accessions are taken from the preprint Data availability section. Confirm that each GSE is public at submission time.

Full table: [`reproducibility/data_accessions.tsv`](reproducibility/data_accessions.tsv).

---

## System requirements

| Component | Status | Notes |
|-----------|--------|-------|
| OS | **TODO (author)** | Scripts were developed on macOS (`/Users/jonathan/...`) and Linux (`/mnt/large/...`). State the exact OS/distro tested for reviewers. |
| R | **Partial** | Trajectory scripts: **R 4.2.2** (script comment). SAMMY visualization stack in Methods: **R 4.4.2**. Bulk DESeq2 **1.30.1** implies an older Bioconductor than the R 4.4.2 stack — expect **multiple environments**. |
| Python | **Partial** | QC notebook metadata: **Python 3.7.7** (`xci-venv`). |
| Nextflow | **TODO (author)** | Required for nf-core/sammyseq and nf-core/chipseq; version not recorded in repo. |
| Conda/Mamba | Optional | Recommended once authors export env files. |
| External CLI | STAR 2.7.10a, cutadapt 4.1, Bowtie2 2.2.5, Cell Ranger 6.1.2, Salmon (version TODO), FastQC/MultiQC | Mostly upstream of deposited scripts |
| Hardware | **TODO (author)** | `tradeSeq` script defaults to 24 cores; WebGestaltR uses 20 threads. Laptop may suffice for demo and DESeq2; full scRNA/SAMMY needs more RAM (**TODO measure**). |
| Non-standard hardware | None known beyond multi-core CPU | No GPU requirement identified in scripts |

See [`reproducibility/software_versions.tsv`](reproducibility/software_versions.tsv).

---

## Installation

Exact bit-for-bit recreation is **not** possible from this repository alone (no `renv.lock` / conda lockfile). Use author-exported environments when available. Inventories:

- [`reproducibility/environments/bulk_r_packages.txt`](reproducibility/environments/bulk_r_packages.txt)
- [`reproducibility/environments/scrnaseq_r_packages.txt`](reproducibility/environments/scrnaseq_r_packages.txt)
- [`reproducibility/environments/python_packages.txt`](reproducibility/environments/python_packages.txt)
- [`reproducibility/environments/nextflow_pins.env.example`](reproducibility/environments/nextflow_pins.env.example)

### Demo-only (Python 3, no R)

```bash
git clone https://github.com/tartaglialabIIT/LBROmics.git
cd LBROmics
# Uses stdlib + optional matplotlib if installed
python3 demo/run_demo.py
```

**Installation time (demo):** typically **&lt; 1 minute** to clone; demo dependencies are Python standard library (matplotlib optional).  
**Full analysis stack install time:** **TODO (author measure)** on a normal workstation after providing lockfiles.

### Bulk RNA-seq (R)

1. Install an R version compatible with **DESeq2 1.30.1** (prefer restoring from author `sessionInfo()` / `renv.lock` — **TODO author**).
2. Install packages listed in `reproducibility/environments/bulk_r_packages.txt`.
3. Place GEO count matrices and GRCm38.98 GTF as described in [`reproducibility/execution_order.md`](reproducibility/execution_order.md).

### scRNA-seq

1. Python **3.7.7**-compatible env with Scanpy **1.9.1** (prefer `pip freeze` from `xci-venv` — **TODO author**).
2. Separate R **4.2.2** env with Seurat **4.1.0** and trajectory packages (prefer `sessionInfo()` — **TODO author**).

### 4f-SAMMY-seq / ChIP-seq

1. Install Nextflow (**version TODO**).
2. Pin nf-core/sammyseq and nf-core/chipseq commits (**TODO author**; Methods cite sammyseq `dev`).
3. Use R **4.4.2** for Gviz-based visualization if reproducing Methods figures beyond the two deposited R scripts.

---

## Demo

A small demo validates that deposited bulk DE tables are readable and produces a volcano-style summary plot plus a sanity-check against manuscript DEG counts.

```bash
python3 demo/run_demo.py
# or: make demo
```

| Item | Value |
|------|-------|
| Input | `bulk/deseq2_results_NPC.txt`, `bulk/deseq2_results_mESC.txt` (already in repo) |
| Command | `python3 demo/run_demo.py` |
| Expected output | `demo/output/demo_summary.txt`, `demo/output/volcano_NPC.png` (PNG if matplotlib available; otherwise summary only) |
| Approximate runtime | ~0.1 s (summary only) or ~2 s with matplotlib volcano on a Linux workstation (measured 2026-09-21) |

See [`demo/README.md`](demo/README.md).

---

## Instructions for use

### Expected inputs

| Workflow | Inputs |
|----------|--------|
| Bulk DE | Tab-delimited raw counts (`gene_id` + 12 samples), Ensembl GRCm38.98 GTF |
| Karyoplot | Deposited `deseq2_results_NPC.txt` + `genes_for_karyo.csv` (included) |
| WebGestalt | DESeq2 result tables with `gene_id`, `log2FoldChange`, `padj` |
| scRNA | 10x MTX directories; `good_cells.csv` from QC |
| SAMMY viz | Differential gene-list files / analysis RDS + mm10 GTF + Salmon gene TPM matrix |

### Execution order

See [`reproducibility/execution_order.md`](reproducibility/execution_order.md).

Briefly:

1. **Bulk:** `deseq2_analysis.R` → `webgestalt.R` → `karyoplot_mouse.R`
2. **scRNA:** notebook `1_` → `2_` → `3_` → `4_` → `5_` → `6_` → `7_`
3. **SAMMY:** run nf-core/sammyseq externally → `Barplot_DiffSolGenes.R` / `TPM_DSR.R`

Path overrides (after portability updates): set `LBROMICS_SCRNA_ROOT`, `LBROMICS_SAMMY_ROOT`, or pass CLI arguments where supported — see script headers.

### Own data

- **Bulk:** provide count matrices with the same layout (or edit the sample-name mapping in `deseq2_analysis.R`; do not rely on silent column-order renaming for new designs).
- **scRNA:** point QC/clustering scripts at your 10x directories; re-tune QC thresholds only if scientifically justified.
- **SAMMY viz:** point scripts at your differential-solubility outputs.

### Figure ↔ code map

[`reproducibility/figure_to_code.tsv`](reproducibility/figure_to_code.tsv) — partial; panels marked `TODO` need author confirmation.

---

## Code availability (manuscript-ready text)

> Custom R and Python scripts used to reproduce the bulk RNA-seq, single-cell RNA-seq, and 4f-SAMMY-seq visualization analyses are publicly available at https://github.com/tartaglialabIIT/LBROmics (MIT license). A versioned archive of this repository will be deposited in Zenodo upon acceptance (**TODO: add Zenodo DOI after tagging a release**). Sequencing data are available in GEO under accessions GSE318892, GSE318871, GSE318873, GSE318895, GSE324396, with related public datasets GSE165447 and GSE96107.

---

## Citation

Please cite the preprint/paper and this repository. Machine-readable metadata: [`CITATION.cff`](CITATION.cff).

Fiorentino J, Perotti I, Ruiz Blanes N, Rosti V, Sigala I, Nikolakaki E, Colantoni A, D’Elia A, Massari R, Scavizzi F, Raspa M, Ascolani M, Humphreys NE, Giannakouros T, Guttman M, Lanzuolo C, Tartaglia GG, Cerase A. LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization. bioRxiv doi:10.64898/2026.03.30.714681.

---

## Smoke tests

```bash
make check
# or
bash scripts/smoke_test.sh
```

---

## Abstract (from manuscript materials)

The nuclear lamina plays a central role in genome organization, yet how specific lamina-associated proteins regulate chromosome architecture during development remains unclear. Here, we show that the nucleoplasmic domains of the Lamin B Receptor (LBR) are essential for X-chromosome localization at the nuclear periphery and chromatin architecture during neural differentiation. Using genetic dissection of LBR function, combined with genome-wide chromatin solubility profiling and transcriptional analyses, we demonstrate that loss of LBR N-terminal domains impairs proper cell differentiation and X chromosome inactivation (XCI), selectively disrupting chromatin structure in neural progenitors but not in pluripotent cells.

For the peer-reviewed wording of X-chromosome solubility phenotypes, refer to the preprint/journal version (the repository README previously contained a simplified abstract that may differ in solubility direction wording — see audit).
