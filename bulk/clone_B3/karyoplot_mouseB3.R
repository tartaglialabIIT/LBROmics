library(biomaRt)
library(regioneR)
library(karyoploteR)

# Expect DE results and genes_for_karyo_B3.csv in the working directory
# (or setwd() to your DiffExp output folder)

genes.for.karyo <- read.csv("genes_for_karyo_B3.csv")
genes.for.karyo <- genes.for.karyo$gene_name

DE.res <- read.csv("deseq2_results_NPC.csv")
DE.res.karyo <- DE.res[DE.res$gene_name %in% genes.for.karyo, ]
DE.res.karyo.upWT <- DE.res.karyo$gene_name[DE.res.karyo$log2FoldChange < 0]
DE.res.karyo.upMut <- DE.res.karyo$gene_name[DE.res.karyo$log2FoldChange > 0]

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

pdf("karyoplot_Xchr_clone_B3.pdf", width = 14, height = 9)
kp <- plotKaryotype(genome = "mm10", chromosomes = c("chrX"), plot.type = 2)
kpPlotMarkers(kp, data = karyo.upWT, labels = karyo.upWT$external_gene_name,
              text.orientation = "vertical", adjust.label.position = TRUE,
              data.panel = 1, label.color = "black", max.iter = 1000)
kpPlotMarkers(kp, data = karyo.upMut, labels = karyo.upMut$external_gene_name,
              text.orientation = "vertical", adjust.label.position = TRUE,
              data.panel = 2, label.color = "black", max.iter = 1000)
dev.off()
