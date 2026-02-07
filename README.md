# LBROmics

*LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization*

## Abstract

The nuclear lamina plays a central role in genome organization, yet how specific lamina-associated proteins regulate chromosome architecture during development remains unclear. Here, we show that the nucleoplasmic domains of the Lamin B Receptor (LBR) are essential for X-chromosome localization at the nuclear periphery and chromatin architecture during neural differentiation. Using genetic dissection of LBR function, combined with genome-wide chromatin solubility profiling and transcriptional analyses, we demonstrate that loss of LBR N-terminal domains impairs proper cell differentiation and X chromosome inactivation (XCI), selectively disrupting chromatin structure in neural progenitors but not in pluripotent cells.

Strikingly, these effects are disproportionately concentrated on the inactive X chromosome, which undergoes a pronounced shift toward an insoluble chromatin state, revealing a decoupling between chromatin solubility and steady-state gene expression. Our findings establish the nucleoplasmic function of LBR as a key determinant of X-chromosome functionality and identify chromatin solubility and accessibility as a previously underappreciated dimension of genome regulation by the nuclear lamina in XCI. Finally, our work provides definitive genetic evidence that LBR’s nuclear architectural functions are molecularly separable from its metabolic sterol reductase activity, which is preserved in our model, and are critically necessary for XCI in differentiating mouse female XX ESCs models.

## Repository description

This repository contains the analysis code used in the study “LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization.”

It provides scripts and workflows to reproduce the main computational analyses of the manuscript, including bulk RNA-seq, single-cell RNA-seq, and chromatin solubility (4f-SAMMY-seq) analyses.

## Repository structure

**bulkRNAseq/**

Scripts for reproducing bulk RNA-seq analyses, including:
	•	differential expression analysis (DESeq2)
	•	PCA and sample distance quality controls
	•	heatmaps for key gene sets (escapees, X-linked genes, markers)
	•	annotated karyoplot of X chromosome
	•	GSEA and ORA using webgestaltR

These scripts use raw count matrices deposited on GEO.

**scRNAseq/**

Scripts for single-cell RNA-seq analysis, including:
	•	preprocessing and clustering
	•	differential expression analyses
	•	marker gene visualization
	•	trajectory inference and tradeseq analysis

**4fSAMMYseq/**

Analysis workflows for chromatin solubility profiling (4f-SAMMY-seq), including:
	•	differential solubility analysis
	•	chromatin compartment analysis
	•	comparison with ChIP-seq data

## Data availability

Raw and processed sequencing data are available on GEO (accession numbers provided in the manuscript).

## Citation

If you use this code or build upon these analyses, please cite:

Fiorentino J†, Perotti I†, Ruiz Blanes N, Rosti V, Gammon L, Sigala I, Nikolakaki E, Colantoni A, D’Elia A, Massari R, Scavizzi F, Raspa M, Ascolani M, Humphreys NE, Giannakouros T, Guttman M, Lanzuolo C, Tartaglia GG, Cerase A.
LBR nucleoplasmic domains regulate X-chromosome solubility and nuclear organization.
