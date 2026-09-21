# LBROmics

Analysis code for:

> Fiorentino et al. *LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization.*  
> bioRxiv [doi:10.64898/2026.03.30.714681](https://doi.org/10.64898/2026.03.30.714681)

Custom R and Python workflows for bulk RNA-seq, single-cell RNA-seq, and 4f-SAMMY-seq analyses described in the manuscript.

**License:** MIT · **Repository:** https://github.com/tartaglialabIIT/LBROmics

---

## Contents

| Directory | Description |
|-----------|-------------|
| `bulk/` | Preprocessing, DESeq2, enrichment, karyoplot, chromosomal enrichment, liver DE, clone B3 helpers |
| `scRNAseq/` | Cell Ranger helper, Scanpy QC, Seurat clustering/markers/trajectories, tradeSeq, CellRank/scvelo, in silico bulk integration |
| `4fSAMMYseq/` | nf-core samplesheet/pins usage, DiffSol barplots, TPM–DSR, compartments, ChIP comparison |
| `demo/` | Small demo on deposited DE tables |
| `reproducibility/` | GEO accessions, software versions, environments, run order |
| `CITATION.cff` | Citation metadata |

---

## Data availability

| Analysis | GEO | Code |
|----------|-----|------|
| Bulk RNA-seq (clone A8) | [GSE318892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318892) | `bulk/` |
| Bulk RNA-seq (clone B3) | [GSE318895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318895) | `bulk/clone_B3/` |
| scRNA-seq | [GSE318871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318871) | `scRNAseq/` |
| 4f-SAMMY-seq | [GSE318873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318873) | `4fSAMMYseq/` |
| Liver bulk RNA-seq | [GSE324396](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324396), [GSE165447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165447) | `bulk/liver_young_kumar/` |
| Public ChIP-seq (Bonev et al.) | [GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107) | `4fSAMMYseq/` |

See also `reproducibility/data_accessions.tsv`.

---

## Requirements

Key versions are listed in `reproducibility/software_versions.tsv`. In brief:

- Python **3.7.7** with packages in `reproducibility/environments/python_xci-venv/requirements.txt`
- R **4.2.2** (scRNA-seq trajectories) and R **4.4.2** (some SAMMY visualization)
- Nextflow **24.10.4**; nf-core/sammyseq **`dev` @ `fa6f6ffeb3`** (`reproducibility/environments/nextflow_pins.env`)
- Genome annotation: Ensembl GRCm38.98 / mm10

---

## Quick start

```bash
git clone https://github.com/tartaglialabIIT/LBROmics.git
cd LBROmics
python3 demo/run_demo.py
```

Full workflows: [`reproducibility/HOW_TO_RUN.md`](reproducibility/HOW_TO_RUN.md).

---

## Code availability

Custom R and Python scripts used for the bulk RNA-seq, single-cell RNA-seq and 4f-SAMMY-seq analyses are available at https://github.com/tartaglialabIIT/LBROmics (MIT license). Sequencing data are available in GEO (accessions above).

---

## Citation

See [`CITATION.cff`](CITATION.cff).
