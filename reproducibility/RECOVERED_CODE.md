# Recovered analysis scripts

Many scripts retain historical absolute paths (`/Users/...`, `/mnt/large/...`, `/mnt/storage/...`). Adapt paths or set working directories before running.

## Bulk (`bulk/`)

| Path | Role |
|------|------|
| `preprocessing/` | cutadapt, FastQC, rRNA depletion, STAR alignment helpers |
| `liver_young_kumar/` | Liver DE comparing Young et al. / Kumar samples |
| `chromosomal_enrichment/` | Chromosome enrichment notebook (Python) |
| `clone_B3/` | Clone B3 DE context + karyoplot companions |

## scRNA-seq (`scRNAseq/`)

| Path | Role |
|------|------|
| `cellranger_count_script.sh` | Cell Ranger count commands |
| `RNA_velocity/` | kb_python launch + scvelo/CellRank notebook |
| `in_silico_bulk_DE/` | Pseudobulk DE from single-cell counts |
| `Integrated_PCA/` | Integration of in silico bulk with bulk RNA-seq |
| `good_cells.csv` | QC-pass cell list for Seurat clustering (see note below) |

`good_cells.numbers` is a legacy Apple Numbers document. Prefer `good_cells.csv` (3560 cells). Confirm against your original Numbers export if needed.

## 4f-SAMMY-seq (`4fSAMMYseq/`)

| Path | Role |
|------|------|
| `pipeline_launch/` | Historical nf-core/sammyseq samplesheets |
| `compartments_analysis.R` | Compartment visualization |
| `downstream_viz/` | ChIP comparison, Gviz tracks, DiffSol extras |
| `ChIP_peaks_vs_random/` | Peak enrichment vs matched random regions |

Launch pins: `reproducibility/environments/nextflow_pins.env` and `reproducibility/environments/nextflow/`.
