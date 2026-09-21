# LBROmics

Analysis code accompanying:

> Fiorentino et al. *LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization.*  
> bioRxiv [doi:10.64898/2026.03.30.714681](https://doi.org/10.64898/2026.03.30.714681)

This repository provides custom R and Python workflows for bulk RNA-seq, single-cell RNA-seq, and downstream 4f-SAMMY-seq visualization.

**License:** MIT  
**Repository:** https://github.com/tartaglialabIIT/LBROmics  
**Reproducibility materials:** [`reproducibility/`](reproducibility/)

---

## Repository structure

| Path | Contents |
|------|----------|
| `bulk/` | DESeq2 analysis, QC/heatmaps, WebGestaltR enrichment, X-chromosome karyoplot; deposited DE result tables |
| `scRNAseq/` | Scanpy QC notebook, Seurat clustering/markers, destiny/slingshot/condiments, tradeSeq, gProfiler |
| `4fSAMMYseq/` | Visualization of differential solubility gene counts and TPM in differentially soluble regions |
| `demo/` | Small runnable demo using deposited DE tables |
| `reproducibility/` | Data accessions, software inventory, figure map, execution order, environment notes |
| `scripts/` | Smoke tests and helpers |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | MIT |

---

## Data availability

| Analysis | Biological system | GEO accession | Repo directory |
|----------|-------------------|---------------|----------------|
| Bulk RNA-seq (clone A8) | Female XX mESC / day-5 NPC | [GSE318892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318892) | `bulk/` |
| Bulk RNA-seq (clone B3) | Female XX mESC / day-5 NPC | [GSE318895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318895) | see manuscript |
| scRNA-seq | WT vs Lbr NT-KO day-5 NPC | [GSE318871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318871) | `scRNAseq/` |
| 4f-SAMMY-seq | ESC / NPC clone B3 | [GSE318873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318873) | `4fSAMMYseq/` |
| Bulk RNA-seq liver (new) | Male/female liver | [GSE324396](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324396) | see manuscript |
| Bulk RNA-seq liver (prior) | Young et al. subset | [GSE165447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165447) | see manuscript |
| ChIP-seq (public) | Bonev et al. NPCs | [GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107) | see manuscript |

Full table: [`reproducibility/data_accessions.tsv`](reproducibility/data_accessions.tsv).

---

## System requirements

Versions used in the study are summarized in [`reproducibility/software_versions.tsv`](reproducibility/software_versions.tsv) and in the manuscript Methods. In brief:

- **R:** multiple analysis environments were used (trajectory scripts note R 4.2.2; some SAMMY visualization used R 4.4.2)
- **Python:** 3.7.7 for Scanpy QC (`scRNAseq/1_QualityControl.ipynb` metadata); Scanpy 1.9.1
- **Key R packages:** DESeq2 1.30.1, Seurat 4.1.0, WebGestaltR 0.4.5, and others listed in `reproducibility/environments/`
- **Reference:** Ensembl GRCm38.98 / mm10
- **Hardware:** multi-core CPU recommended for tradeSeq / WebGestaltR (`nThreads` / parallel workers are set in those scripts)

Package inventories (not lockfiles):

- [`reproducibility/environments/bulk_r_packages.txt`](reproducibility/environments/bulk_r_packages.txt)
- [`reproducibility/environments/scrnaseq_r_packages.txt`](reproducibility/environments/scrnaseq_r_packages.txt)
- [`reproducibility/environments/python_packages.txt`](reproducibility/environments/python_packages.txt)
- [`reproducibility/environments/nextflow_pins.env.example`](reproducibility/environments/nextflow_pins.env.example)

---

## Installation

```bash
git clone https://github.com/tartaglialabIIT/LBROmics.git
cd LBROmics
```

For the demo only, Python 3 is sufficient (matplotlib optional).

For full analyses, install R/Python packages from the inventories above, preferably matching the versions reported in the manuscript Methods / `software_versions.tsv`.

---

## Demo

```bash
python3 demo/run_demo.py
# or: make demo
```

| Item | Value |
|------|-------|
| Input | `bulk/deseq2_results_NPC.txt`, `bulk/deseq2_results_mESC.txt` |
| Expected output | `demo/output/demo_summary.txt`, optional `demo/output/volcano_NPC.png` |
| Runtime | ~0.1–2 s on a typical workstation |

See [`demo/README.md`](demo/README.md).

---

## Instructions for use

Step-by-step execution order, inputs, and path environment variables:

- [`reproducibility/execution_order.md`](reproducibility/execution_order.md)
- Figure-oriented map: [`reproducibility/figure_to_code.tsv`](reproducibility/figure_to_code.tsv)

**Bulk (summary):** place GEO count matrices and GRCm38.98 GTF as described in `execution_order.md`, then run `deseq2_analysis.R`, `webgestalt.R`, and `karyoplot_mouse.R`.

**scRNA-seq (summary):** numbered scripts `1_`…`7_` in `scRNAseq/`. Set `LBROMICS_SCRNA_ROOT` (and related env vars) as described in the script headers.

**4f-SAMMY-seq visualization (summary):** set `LBROMICS_SAMMY_ROOT` or pass `--folder` / `--rdata` / `--gtf` / `--tpm` arguments to the scripts in `4fSAMMYseq/`.

---

## Code availability (manuscript text)

> Custom R and Python scripts used for the bulk RNA-seq, single-cell RNA-seq, and 4f-SAMMY-seq visualization analyses are publicly available at https://github.com/tartaglialabIIT/LBROmics (MIT license). Sequencing data are available in GEO under accessions GSE318892, GSE318871, GSE318873, GSE318895, GSE324396, with related public datasets GSE165447 and GSE96107.

---

## Citation

See [`CITATION.cff`](CITATION.cff).

Fiorentino J, Perotti I, Ruiz Blanes N, Rosti V, Sigala I, Nikolakaki E, Colantoni A, D’Elia A, Massari R, Scavizzi F, Raspa M, Ascolani M, Humphreys NE, Giannakouros T, Guttman M, Lanzuolo C, Tartaglia GG, Cerase A. LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization. bioRxiv doi:10.64898/2026.03.30.714681.

---

## Smoke tests

```bash
make check
```
