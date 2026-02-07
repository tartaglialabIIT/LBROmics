library(biomaRt)
library(regioneR)

wd <- "/mnt/large/jfiorentino/karyoplot"
setwd(wd)

# Read the genes to be highlighted
genes.for.karyo <- read.csv("genes_for_karyo.csv")
genes.for.karyo <- genes.for.karyo$gene_name

# Load the results from differential expression
DE.res <- read.csv("./deseq2_results_NPC.txt",sep="\t")
DE.res.karyo <- DE.res[DE.res$gene_name %in% genes.for.karyo,]
DE.res.karyo.upWT <- DE.res.karyo[DE.res.karyo$log2FoldChange<0,]
DE.res.karyo.upWT <- DE.res.karyo.upWT$gene_name
DE.res.karyo.upMut <- DE.res.karyo[DE.res.karyo$log2FoldChange>0,]
DE.res.karyo.upMut <- DE.res.karyo.upMut$gene_name

ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

karyo.upWT <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
               filters = 'external_gene_name', values =DE.res.karyo.upWT, mart = ensembl))
seqlevelsStyle(karyo.upWT) <- "UCSC"

karyo.upMut <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
                              filters = 'external_gene_name', values =DE.res.karyo.upMut, mart = ensembl))
seqlevelsStyle(karyo.upMut) <- "UCSC"

color1 <- '#F8766D'
color2 <- '#00BFC4'

library(karyoploteR) 
pdf("mm10_DE_Xchr_sig_genes.pdf",width=14,height=9)


kp <- plotKaryotype(genome="mm10",chromosomes = c("chrX"),plot.type=2)
kpPlotMarkers(kp, data=karyo.upWT, labels=karyo.upWT$external_gene_name, text.orientation = "vertical",
              #r1=1.1, cex=0.8, 
              adjust.label.position = T,data.panel=1,label.color='black', max.iter = 1000)
kpPlotMarkers(kp, data=karyo.upMut, labels=karyo.upMut$external_gene_name, text.orientation = "vertical",
              #r1=1.1, cex=0.8, 
              adjust.label.position = T,data.panel=2,label.color='black',max.iter = 1000)

dev.off()
