# Code & reproducibility audit

**Repository:** https://github.com/tartaglialabIIT/LBROmics  
**Manuscript / preprint:** Fiorentino et al., *LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization*  
**bioRxiv DOI:** [10.64898/2026.03.30.714681](https://doi.org/10.64898/2026.03.30.714681)  
**Audit date:** 2026-09-21  
**Commit audited:** `e032ef9` (pre-improvement `main`)

This audit separates **what is already in the repository** from **what can be recovered from the preprint methods** vs **what still requires author input**. Version numbers taken from the preprint are cited explicitly; they are **not** package lockfiles and do **not** alone recreate identical environments.

---

## Summary verdict

| Area | Can exact environment be recreated from the repo alone? |
|------|---------------------------------------------------------|
| bulk RNA-seq (DESeq2 / WebGestalt / karyoplot) | **No** — no lockfile; some versions in preprint; hard-coded absolute paths |
| scRNA-seq (Scanpy QC → Seurat → destiny/slingshot/tradeSeq) | **No** — no lockfile; many versions in preprint; absolute paths; intermediate objects missing |
| 4f-SAMMY-seq (nf-core/sammyseq + viz) | **No** — pipeline configs/commits missing; only post-hoc visualization scripts present |
| ChIP-seq reprocessing (nf-core/chipseq) | **No** — no scripts or pinned pipeline version in repo |
| Liver bulk RNA-seq / chromosomal enrichment Python | **No** — described in manuscript, **not present** in repo |

---

## Reproducibility / journal code-sharing checklist audit table

| Requirement | Current status | Evidence found in repo | What should be added | Can this be inferred safely from the repo? | User input needed? |
|-------------|----------------|------------------------|----------------------|--------------------------------------------|--------------------|
| Source code present | PASS | R scripts + Scanpy notebook under `bulk/`, `scRNAseq/`, `4fSAMMYseq/` | Missing manuscript analyses (see gaps) | yes | partial (to deposit missing analysis code) |
| Small demo / test dataset | MISSING (pre-improvement) | No `demo/`; no toy inputs; `good_cells.numbers` is Apple Numbers binary, not usable CSV | Runnable demo + tiny inputs + expected output + runtime | no | if full-workflow demo needed |
| System requirements: OS | MISSING | Paths imply macOS (`/Users/jonathan/...`) and Linux (`/mnt/large/...`); no tested-OS statement | Explicit tested OS/distro | no | **yes** |
| System requirements: R version | PARTIAL | Comment in `scRNAseq/5_trajectory_inference_and_condiments.R`: “runs in R4.2.2”; preprint Methods: R 4.4.2 used for some SAMMY visualization | Single tested matrix per workflow; reconcile multi-R usage | yes (partial from script + preprint) | **yes** (confirm which R for each branch) |
| System requirements: Python version | PARTIAL | Notebook metadata: Python **3.7.7**, kernel `xci-venv` (`scRNAseq/1_QualityControl.ipynb`) | Confirm whether 3.7.7 was used for all Python steps (scipy/statsmodels/scvelo/cellrank) | yes (notebook only) | **yes** |
| System requirements: Nextflow | MISSING in repo | Preprint: nf-core/sammyseq (`/dev/`) and nf-core/chipseq used; **no** Nextflow configs in repo | Nextflow version, pipeline release/commit, profile, configs | no (preprint only, unpinned) | **yes** |
| Dependencies with version numbers | PARTIAL | Preprint Methods list many versions; repo has **no** `renv.lock`, `environment.yml`, `requirements.txt`, containers, or `sessionInfo()` | Author-exported lockfiles / `sessionInfo()` / conda exports | versions: yes from preprint; lockfile: **no** | **yes** |
| Tested versions / platforms | MISSING | No “tested on” section | Document machines/OS/R/Python actually used | no | **yes** |
| Non-standard hardware | PARTIAL | `6_run_tradeseq.R` sets `BPPARAM$workers <- 24`; `webgestalt.R` uses `nThreads=20` | State CPU/RAM expectations; whether HPC required | partial | **yes** (RAM/walltime) |
| Installation guide | MISSING | README has abstract/structure only | Exact install steps per branch | no | **yes** (until lockfiles exist) |
| Typical installation time | MISSING | — | Measured install time on a desktop | no | **yes** |
| Demo: commands, expected output, runtime | MISSING | — | Demo script + expected artefacts | no | preferred |
| Instructions for use / own data | PARTIAL | `deseq2_analysis.R` documents expected count-table format; other scripts assume local absolute paths | Path config, input formats, execution order per workflow | partial | **yes** |
| Reproduce manuscript analyses | PARTIAL | Numbered scRNA scripts suggest order; bulk/SAMMY incomplete vs Methods | End-to-end commands + figure map + GEO file mapping | partial | **yes** |
| Figure ↔ script mapping | MISSING | Filenames (e.g. `markers_figs6.pdf`) give weak hints only | Curated `figure_to_code.tsv` with author review | no (safe mapping incomplete) | **yes** |
| GEO accession mapping | PARTIAL in repo / PASS in preprint | README: “accession numbers provided in the manuscript”; preprint Data availability lists GSE accessions | Table in README / `data_accessions.tsv` | yes from preprint | confirm public status |
| Exact pipeline versions / commits | MISSING | Preprint cites nf-core/sammyseq **dev** and compartments branch URL; no commit hash; chipseq version absent | Pin release or `git rev-parse` for sammyseq & chipseq | no | **yes** |
| Output locations documented | PARTIAL | Scripts write under hard-coded local dirs (`plots_JULY_2023/`, `/mnt/large/...`) | Relative `results/` layout documented | partial | no (docs can fix) |
| License | PASS | MIT `LICENSE` (Copyright 2026 tartaglialabIIT) | Keep; mention in README | yes | no |
| Citation / manuscript link | PARTIAL | Authors in README; no DOI/preprint link; no `CITATION.cff` | `CITATION.cff` + preprint DOI + Zenodo TODO | yes | ORCID optional |
| Repository URL | PASS | GitHub public | Also in manuscript code availability | yes | no |
| Version / release | MISSING | No git tags / GitHub release / Zenodo DOI | Tag release + archive on Zenodo | no | **yes** |
| Container images / digests | MISSING | None | Docker/Singularity if used with nf-core | no | **yes** |
| SLURM / HPC config | MISSING | Absolute `/mnt/large/` paths suggest shared filesystem | Document executor if any | no | **yes** |

---

## Environment artefact inventory (repo search)

| Artefact | Present? |
|----------|----------|
| `environment.yml` / `.yaml` | **No** |
| `requirements.txt` / `pyproject.toml` / `poetry.lock` | **No** |
| `pip freeze` output | **No** |
| `renv.lock` / packrat | **No** |
| R `DESCRIPTION` | **No** |
| `sessionInfo()` dump | **No** |
| Bioconductor version file | **No** |
| Conda explicit export | **No** |
| Nextflow config / profiles | **No** |
| nf-core pipeline pin (version/commit) | **No** |
| Docker / Singularity / Apptainer | **No** |
| Shell driver scripts with execution order | **No** |
| `.gitignore` | **No** (pre-improvement) |

**Versions found outside lockfiles (safe to record, not sufficient to recreate):**

| Software | Version | Source |
|----------|---------|--------|
| Python | 3.7.7 | `scRNAseq/1_QualityControl.ipynb` metadata |
| R | 4.2.2 | comment in `scRNAseq/5_trajectory_inference_and_condiments.R` |
| R | 4.4.2 | preprint Methods (SAMMY visualization with Gviz/rtracklayer) |
| DESeq2 | 1.30.1 | preprint Methods |
| WebGestaltR | 0.4.5 | preprint Methods |
| pheatmap | 1.0.12 | preprint Methods |
| karyoploteR | 1.20.3 | preprint Methods |
| cutadapt | 4.1 | preprint Methods |
| Bowtie2 | 2.2.5 | preprint Methods |
| STAR | 2.7.10a | preprint Methods |
| scipy | 1.7.2 | preprint Methods |
| statsmodels | 0.13.1 | preprint Methods |
| Cell Ranger | 6.1.2 (GRCm38) | preprint Methods |
| Scanpy | 1.9.1 | preprint Methods |
| Seurat | 4.1.0 | preprint Methods |
| limma | 3.46.0 | preprint Methods |
| destiny | 3.4.0 | preprint Methods |
| slingshot | 2.8.0 | preprint Methods |
| condiments | 1.6.0 | preprint Methods |
| tradeSeq | 1.8.0 | preprint Methods |
| scry | 1.6.0 | preprint Methods |
| gprofiler2 | 0.2.1 (trajectory DEGs); 0.2.4 (DSR GO) | preprint Methods |
| kb_python | 0.27.3 | preprint Methods |
| scvelo | 0.2.5 | preprint Methods |
| cellrank | 1.5.1 | preprint Methods |
| Gviz | 1.50.0 | preprint Methods |
| rtracklayer | 1.66.0 | preprint Methods |
| FastQC | 0.12.1 | preprint Methods |
| MultiQC | 1.31 | preprint Methods |
| ImageJ | 1.54f | preprint Methods (not omics pipeline) |
| Genome annotation | Ensembl GRCm38.98 / mm10 | scripts + preprint |
| nf-core/sammyseq | **dev** (unpinned) + `compartments_subworkflow` branch | preprint Methods |
| nf-core/chipseq | **unspecified** | preprint Methods |

**Important:** DESeq2 1.30.1 / limma 3.46.0 / destiny 3.4.0 / slingshot 2.8.0 imply **older Bioconductor releases** than Gviz 1.50.0 / rtracklayer 1.66.0 (R 4.4.x). Multiple computational environments were almost certainly used. A single lockfile would be misleading without author confirmation.

---

## Per-workflow notes

### bulk/

| Item | Status |
|------|--------|
| Scripts | `deseq2_analysis.R`, `webgestalt.R`, `karyoplot_mouse.R` |
| Deposited intermediates | `deseq2_results_mESC.txt`, `deseq2_results_NPC.txt`, `genes_for_karyo.csv` |
| Raw inputs in repo | **No** (`raw_counts_*.txt`, GTF) — expected from GEO |
| Paths | `deseq2_analysis.R` uses relative paths; `webgestalt.R` / `karyoplot_mouse.R` use `/mnt/large/jfiorentino/...` |
| Missing vs Methods | Liver DE; chromosomal enrichment (Python scipy/statsmodels); STAR/cutadapt/Bowtie2 preprocessing; clone B3 bulk |

### scRNAseq/

| Item | Status |
|------|--------|
| Order | Numbered `1_`…`7_` implies workflow order |
| QC | Scanpy notebook; placeholder `main_folder='/pathtoyourfolder/'` |
| Clustering onward | Hard-coded `/Users/jonathan/Desktop/...` and `/mnt/large/...` |
| `good_cells.numbers` | Apple Numbers document — **not** the `good_cells.csv` expected by `2_clustering.R` |
| Missing vs Methods | CellRanger outputs; CellRank/scvelo/kb_python; in silico bulk + limma integration; RDS intermediates |

### 4fSAMMYseq/

| Item | Status |
|------|--------|
| Present | Visualization only: `Barplot_DiffSolGenes.R`, `TPM_DSR.R` |
| README claim | Differential solubility, compartments, ChIP comparison |
| Reality | Core nf-core/sammyseq run, compartment analysis, ChIP reprocessing **not in repo** |
| Paths | macOS OneDrive absolute paths |

### ChIP-seq

Described in preprint (public GSE96107 + nf-core/chipseq). **No analysis code in this repository.**

---

## Potential scientific issues requiring author review

These were **not** silently changed:

1. **`scRNAseq/3_cluster_markers.R` (approx. lines 29–33):** mutant counts are normalized by `sum.WT` rather than `sum.mutant`:
   ```r
   df[df$cond=="mutant","count"] <- df[df$cond=="mutant","count"]/sum.WT
   ```
   Likely a copy-paste error for the barplot of condition frequencies. Confirm whether published figures used this code or a corrected version.

2. **`scRNAseq/7_gProfiler_lineage2.R` (tail):** references undefined variables `dir.id`, `i`, `my.terms` (looks like leftover from `myfunctions.R`). The early gProfiler block is self-contained; the trailing block will error if run as-is.

3. **Abstract wording mismatch:** repository README abstract says inactive X shifts toward an **insoluble** chromatin state; bioRxiv abstract emphasizes a shift toward a **more soluble** chromatin state (with later discussion of chromosome-specific solubility changes). Align README with the peer-reviewed text.

4. **README vs content for `4fSAMMYseq/`:** directory does not contain differential solubility / compartment / ChIP workflows claimed in the README structure section.

5. **`scRNAseq/good_cells.numbers`:** binary Numbers file cannot feed `2_clustering.R`, which expects `./files/good_cells.csv`. Deposit a CSV export (or regenerate from the QC notebook).

6. **Sample renaming by column order** in `deseq2_analysis.R`: if GEO column order differs from the assumed 12-sample convention, genotypes will be mis-assigned. Prefer an explicit sample-name map when releasing raw counts.

---

## Information I still need from the authors

### Critical for ticking journal code-sharing boxes

1. **Original environment exports** (do not invent):
   - Per machine / workflow: `sessionInfo()` and `BiocManager::version()`
   - `installed.packages()[, c("Package","Version")]` for each R library used
   - Python: `python --version`, `pip freeze` (or `conda env export --from-history` + `conda list --explicit`) for `xci-venv` and any other envs (scvelo/cellrank)
2. **OS / hardware actually tested** (e.g. macOS version, Ubuntu version, CPU cores, RAM).
3. **Installation and demo runtimes** measured on a normal workstation.
4. **Nextflow + nf-core pins:**
   - `nextflow -version`
   - For the checked-out sammyseq and chipseq pipelines: `git rev-parse HEAD`, tag/release if any, and the profile/config used (`-profile docker/singularity/conda`, custom `nextflow.config`)
   - Confirm whether `https://nf-co.re/sammyseq/dev/` and `compartments_subworkflow` are the exact trees used
5. **Container images/digests** if Docker/Singularity were used.
6. **GEO file-level mapping:** which processed files under each GSE correspond to `raw_counts_NPC.txt`, `raw_counts_mESC.txt`, 10x matrices `s1_DR1` / `s2_A8`, Salmon TPM matrix, and SAMMY differential RDS outputs.
7. **Figure-to-script confirmation** for panels that are not obvious from filenames (especially Fig. 1–5 and key supplements).
8. **Missing analysis code** to deposit if reviewers must reproduce Methods sections that are currently absent (liver DE, chromosomal enrichment Python, CellRank/scvelo, ChIP-seq signal aggregation, nf-core launch commands, compartment plots).
9. **Export `good_cells.csv`** (or regenerate) and clarify relationship to `good_cells.numbers`.
10. **Zenodo archive** of a tagged release (DOI) after this branch is reviewed.
11. Confirm whether **Luke Gammon** should appear in formal citation metadata (in README authors list; not in the bioRxiv author list fetched for this audit).

### Suggested recovery commands (run on original analysis machines)

```bash
# Identity of the machine
uname -a
lscpu | head
free -h

# R (in each distinct R library / conda env used)
R --version
Rscript -e 'sessionInfo()'
Rscript -e 'if (requireNamespace("BiocManager", quietly=TRUE)) BiocManager::version()'
Rscript -e 'write.csv(installed.packages()[,c("Package","Version")], "installed_R_packages.csv")'

# Python QC / velocity envs
python --version
pip freeze > pip_freeze.txt
# if conda was used:
conda env export --from-history > environment.from-history.yml
conda list --explicit > conda_explicit.txt

# Nextflow / nf-core
nextflow -version
# inside each pipeline checkout:
git -C /path/to/sammyseq rev-parse HEAD
git -C /path/to/sammyseq describe --tags --always
git -C /path/to/chipseq rev-parse HEAD
```

---

## Gaps relative to manuscript Methods (code not in repo)

| Analysis described in preprint | In repo? |
|--------------------------------|----------|
| Liver bulk DE (GSE165447 + GSE324396) | No |
| cutadapt / Bowtie2-rRNA / STAR preprocessing | No |
| Chromosomal enrichment (scipy / statsmodels) | No |
| Scanpy QC | Yes (notebook) |
| Seurat clustering / markers / trajectories / tradeSeq / gProfiler | Partial (scripts; intermediates missing) |
| CellRank / scvelo / kb_python | No |
| In silico bulk + limma integration | No |
| ChIP-seq signal over DE promoters (GSE96107) | No |
| nf-core/sammyseq differential solubility + compartments | No (only downstream barplot/TPM scripts) |
| nf-core/chipseq reprocessing of Bonev et al. | No |
| Salmon TPM quantification for DSR–expression plots | No (script expects precomputed matrix) |

---

## Post-improvement note

Documentation, inventories, path portability, citation metadata, `.gitignore`, and a **lightweight demo** using deposited DE tables were added to improve reviewer usability. They **do not** invent lockfiles or nf-core commits. Items that require measured runtimes, tested OS confirmation, or author lockfiles remain incomplete until authors supply the information above.
