#!/usr/bin/env Rscript
library(tradeSeq)
library(SingleCellExperiment)

# Set LBROMICS_TRADESEQ_DIR to the directory containing curves.RData, counts.Rdata,
# myconditions.Rdata, and filtered_genes_for_tradeseq.txt
wd <- Sys.getenv("LBROMICS_TRADESEQ_DIR", unset = "")
if (!nzchar(wd)) {
  # Fall back to scRNAseq/ in the repository when present
  wd <- Sys.getenv("LBROMICS_SCRNA_ROOT", unset = getwd())
}
setwd(wd)

load("curves.RData")
load("counts.Rdata")
load("myconditions.Rdata")


my.genes <- read.table("filtered_genes_for_tradeseq.txt")$V2
my.genes <- intersect(my.genes,rownames(counts))

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 24 # use 24 cores
sceGAM <-  fitGAM(counts = counts, # counts matrix
                  sds=curves,
                  conditions= my_conditions,
                  parallel=TRUE,
                  nknots= 6,
                  genes = my.genes,
                  BPPARAM=BPPARAM)


# Save the sce object
saveRDS(sceGAM, "tradeSeq_sce.rds")
