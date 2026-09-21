#!/usr/bin/env Rscript
# Karyoplot of selected significant DE genes on mouse chrX (mm10).
#
# Usage (from repository root or bulk/):
#   Rscript bulk/karyoplot_mouse.R
#   Rscript bulk/karyoplot_mouse.R --de bulk/deseq2_results_NPC.txt --genes bulk/genes_for_karyo.csv --outdir bulk
#
# Environment overrides:
#   LBROMICS_BULK_DIR  directory containing DE table and genes_for_karyo.csv

suppressPackageStartupMessages({
  library(biomaRt)
  library(regioneR)
  library(karyoploteR)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  getwd()
}

bulk_dir <- Sys.getenv("LBROMICS_BULK_DIR", unset = "")
if (!nzchar(bulk_dir)) {
  # Prefer script directory (bulk/) when run via Rscript --file=
  bulk_dir <- script_dir()
}

de_path <- get_arg("--de", file.path(bulk_dir, "deseq2_results_NPC.txt"))
genes_path <- get_arg("--genes", file.path(bulk_dir, "genes_for_karyo.csv"))
outdir <- get_arg("--outdir", bulk_dir)

stopifnot(file.exists(de_path), file.exists(genes_path))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

genes.for.karyo <- read.csv(genes_path)
genes.for.karyo <- genes.for.karyo$gene_name

DE.res <- read.csv(de_path, sep = "\t")
DE.res.karyo <- DE.res[DE.res$gene_name %in% genes.for.karyo, ]
DE.res.karyo.upWT <- DE.res.karyo[DE.res.karyo$log2FoldChange < 0, ]
DE.res.karyo.upWT <- DE.res.karyo.upWT$gene_name
DE.res.karyo.upMut <- DE.res.karyo[DE.res.karyo$log2FoldChange > 0, ]
DE.res.karyo.upMut <- DE.res.karyo.upMut$gene_name

ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

karyo.upWT <- toGRanges(getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = "external_gene_name", values = DE.res.karyo.upWT, mart = ensembl
))
seqlevelsStyle(karyo.upWT) <- "UCSC"

karyo.upMut <- toGRanges(getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = "external_gene_name", values = DE.res.karyo.upMut, mart = ensembl
))
seqlevelsStyle(karyo.upMut) <- "UCSC"

color1 <- "#F8766D"
color2 <- "#00BFC4"

out_pdf <- file.path(outdir, "mm10_DE_Xchr_sig_genes.pdf")
pdf(out_pdf, width = 14, height = 9)

kp <- plotKaryotype(genome = "mm10", chromosomes = c("chrX"), plot.type = 2)
kpPlotMarkers(
  kp, data = karyo.upWT, labels = karyo.upWT$external_gene_name,
  text.orientation = "vertical",
  adjust.label.position = TRUE, data.panel = 1, label.color = "black", max.iter = 1000
)
kpPlotMarkers(
  kp, data = karyo.upMut, labels = karyo.upMut$external_gene_name,
  text.orientation = "vertical",
  adjust.label.position = TRUE, data.panel = 2, label.color = "black", max.iter = 1000
)

dev.off()
message("Wrote ", normalizePath(out_pdf))
