# LBROmics

Analysis code accompanying:

> Fiorentino et al. *LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization.*  
> bioRxiv [doi:10.64898/2026.03.30.714681](https://doi.org/10.64898/2026.03.30.714681)

This repository provides custom R and Python workflows for bulk RNA-seq, single-cell RNA-seq, and 4f-SAMMY-seq analyses and visualization.

**License:** MIT  
**Repository:** https://github.com/tartaglialabIIT/LBROmics  
**Reproducibility materials:** [`reproducibility/`](reproducibility/)

---

## Repository structure

| Path | Contents |
|------|----------|
| `bulk/` | DESeq2, enrichment, karyoplot; preprocessing; liver DE; chromosomal enrichment; clone B3 helpers |
| `scRNAseq/` | Scanpy QC, Seurat clustering/markers/trajectories, Cell Ranger/velocity/CellRank, in silico bulk integration |
| `4fSAMMYseq/` | nf-core launch samplesheets, DiffSol/TPM plots, compartments, ChIP comparison |
| `demo/` | Small runnable demo using deposited DE tables |
| `reproducibility/` | Accessions, software inventory, environments, execution order |
| `scripts/` | Smoke tests |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | MIT |

See also [`reproducibility/RECOVERED_CODE.md`](reproducibility/RECOVERED_CODE.md).

---

## Data availability

| Analysis | Biological system | GEO accession | Repo directory |
|----------|-------------------|---------------|----------------|
| Bulk RNA-seq (clone A8) | Female XX mESC / day-5 NPC | [GSE318892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318892) | `bulk/` |
| Bulk RNA-seq (clone B3) | Female XX mESC / day-5 NPC | [GSE318895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318895) | `bulk/clone_B3/` |
| scRNA-seq | WT vs Lbr NT-KO day-5 NPC | [GSE318871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318871) | `scRNAseq/` |
| 4f-SAMMY-seq | ESC / NPC clone B3 | [GSE318873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318873) | `4fSAMMYseq/` |
| Bulk RNA-seq liver (new) | Male/female liver | [GSE324396](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324396) | `bulk/liver_young_kumar/` |
| Bulk RNA-seq liver (prior) | Young et al. subset | [GSE165447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165447) | `bulk/liver_young_kumar/` |
| ChIP-seq (public) | Bonev et al. NPCs | [GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107) | `4fSAMMYseq/downstream_viz/`, `ChIP_peaks_vs_random/` |

Full table: [`reproducibility/data_accessions.tsv`](reproducibility/data_accessions.tsv).

---

## System requirements

Summarized in [`reproducibility/software_versions.tsv`](reproducibility/software_versions.tsv) and manuscript Methods:

- **R:** multiple stacks (trajectory scripts note R 4.2.2; SAMMY visualization used R 4.4.2)
- **Python:** 3.7.7 (`xci-venv`); see `reproducibility/environments/python_xci-venv/`
- **Nextflow:** 24.10.4; nf-core/sammyseq `dev` @ `fa6f6ffeb3` (docker profile) — `reproducibility/environments/nextflow_pins.env`
- **Reference:** Ensembl GRCm38.98 / mm10
- **Hardware:** multi-core CPU recommended for tradeSeq / WebGestaltR / Nextflow

---

## Installation

```bash
git clone https://github.com/tartaglialabIIT/LBROmics.git
cd LBROmics
```

**Demo only:** Python 3 (matplotlib optional).

**Python analyses (approximate):**

```bash
python3.7 -m venv xci-venv
source xci-venv/bin/activate
pip install -r reproducibility/environments/python_xci-venv/requirements.txt
```

**R analyses:** install packages listed under `reproducibility/environments/*_r_packages.txt`, matching manuscript versions where possible. Exporting `sessionInfo()` from the original machines is recommended (`reproducibility/environments/R_sessionInfo_TODO.md`).

**Nextflow / sammyseq:** see `reproducibility/environments/nextflow/`.

---

## Demo

```bash
python3 demo/run_demo.py
# or: make demo
```

| Item | Value |
|------|-------|
| Input | `bulk/deseq2_results_NPC.txt`, `bulk/deseq2_results_mESC.txt` |
| Expected output | `demo/output/demo_summary.txt`, optional `volcano_NPC.png` |
| Runtime | ~0.1–2 s |

---

## Instructions for use

Full order: [`reproducibility/execution_order.md`](reproducibility/execution_order.md).  
Figure map: [`reproducibility/figure_to_code.tsv`](reproducibility/figure_to_code.tsv).

Many recovered scripts still contain historical absolute paths; adapt them or set the documented environment variables (`LBROMICS_SCRNA_ROOT`, `LBROMICS_SAMMY_ROOT`, `LBROMICS_GOOD_CELLS`, …).

---

## Code availability (manuscript text)

> Custom R and Python scripts used for the bulk RNA-seq, single-cell RNA-seq, and 4f-SAMMY-seq analyses are publicly available at https://github.com/tartaglialabIIT/LBROmics (MIT license). Sequencing data are available in GEO under accessions GSE318892, GSE318871, GSE318873, GSE318895, GSE324396, with related public datasets GSE165447 and GSE96107.

---

## Citation

See [`CITATION.cff`](CITATION.cff).

Fiorentino J, Perotti I, Ruiz Blanes N, Rosti V, Sigala I, Nikolakaki E, Colantoni A, D’Elia A, Massari R, Scavizzi F, Raspa M, Ascolani M, Humphreys NE, Giannakouros T, Guttman M, Lanzuolo C, Tartaglia GG, Cerase A. LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization. bioRxiv doi:10.64898/2026.03.30.714681.

---

## Smoke tests

```bash
make check
```
