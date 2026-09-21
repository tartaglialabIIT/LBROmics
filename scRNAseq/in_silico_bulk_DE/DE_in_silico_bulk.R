wd<- "/mnt/large/jfiorentino/Cerase_Data/scRNA_seq/in_silico_bulk_DE/"
setwd(wd)

library("DESeq2")
library("scran")
library("dplyr")
library("BiocParallel")
register(MulticoreParam(12))

DE_WT_mut <-function(sce){

  dds <- convertTo(sce, type="DESeq2")
  dds$batch <- sce$batch
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  print(dds)
  dds <- DESeqDataSet(dds, design = ~ batch)                     
  dds <- DESeq(dds,parallel=TRUE)
  res <- results(dds, contrast=c("batch","WT","mutant"),parallel=TRUE)
  res <- res[order(res$padj, decreasing = F),]
  res <- res[!is.na(res$padj),]
  res <- as.data.frame.matrix(res)
  res
}

# Read raw data and metadata
data <- read.csv("raw_data.txt",sep = '\t',header = TRUE,row.names = 1)
data <- t(data)
metadata <- read.csv("metadata.csv")
ribo.genes <- rownames(data)[grep("Rp", rownames(data))]
mito.genes <- rownames(data)[grep("mt\\.", rownames(data))]

data <- data[! rownames(data) %in% ribo.genes,]
data <- data[! rownames(data) %in% mito.genes,]

data.sce <- SingleCellExperiment(list(counts=data),
                                 colData=DataFrame(batch=metadata$batch))

res.tot <- DE_WT_mut(data.sce)

write.csv(res.tot,"res_tot.csv", row.names = TRUE)


