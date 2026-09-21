#!/usr/bin/env Rscript
library(Seurat)
library(slingshot)
library(tradeSeq)
library(UpSetR)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(fgsea)
source("myfunctions.R")
library(tidyr)
library(msigdbr)
library(knitr)

# Set LBROMICS_SCRNA_ROOT to the Seurat analysis directory.
wd <- Sys.getenv("LBROMICS_SCRNA_ROOT", unset = "")
if (!nzchar(wd)) wd <- getwd()
setwd(wd)

sce <- readRDS("./files_JULY_2023/tradeSeq_sce.rds")

condRes <- read.csv("./files_JULY_2023/tradeSeq/conditionTest_JULY_2023/condRes.csv",row.names = 1)
conditionGenes_lineage2 <- rownames(condRes)[condRes$padj_lineage2 <= 10^-5]
conditionGenes_lineage2 <- conditionGenes_lineage2[!is.na(conditionGenes_lineage2)]

colnames(condRes)

condRes.lin2 <- condRes[conditionGenes_lineage2,c("waldStat_lineage2","padj_lineage2")]
condRes.lin2 <- condRes.lin2[order(condRes.lin2$padj_lineage2),]
write.csv(file = "DE_genes_lineage2.csv",x = condRes.lin2)

#custom.GO(as.vector(conditionGenes_lineage2),"lineage2","condRes",rownames(sce))

gprofiler_results <- gost(query = as.vector(conditionGenes_lineage2), 
                          organism = "mmusculus", ordered_query = FALSE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = TRUE, 
                          user_threshold = 0.05, correction_method = "g_SCS", 
                          domain_scope = "annotated", custom_bg = rownames(sce), 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE)

mydf <- as.data.frame(gprofiler_results$result)
mydf <- as.data.frame(apply(mydf,2,as.character))
mydf$p_value <- as.numeric(as.character(mydf$p_value))
mydf <- mydf[order(mydf$p_value,decreasing = F),]
write.csv(x=mydf,file="gProfiler_linage2.csv")

mydf.GO <- mydf[mydf$source %in% c("GO:BP","GO:MF", "GO:CC"),]

pdf(file = "ManhattanPlot_Lineage2.pdf", width = 8,height = 6)
# Plot for highlighting specific terms (should be run after gostplot 
# (with interactive=FALSE))
p <- gostplot(gprofiler_results, capped = TRUE, interactive = F)
pp <- publish_gostplot(p,highlight_terms = c(mydf.GO[1:5,"term_id"],"GO:0051276"),
                       width = NA, height = NA, filename = NULL )
pp
dev.off()


# Plot results
p <- gostplot(gprofiler_results, capped = TRUE, interactive = T)
p
my.Dir <- dir.id
dir.create(file.path("./plots_JULY_2023/", my.Dir), showWarnings = FALSE)
full.dir <- paste("./plots_JULY_2023",my.Dir,sep = '/')
file.name <- paste("gprof",i,sep = "_")
file.name <- paste(file.name,"html",sep=".")
# Save the interactive plot to a html file
saveWidget(p, file=paste(full.dir,file.name,sep = "/"))
mydf <- as.data.frame(gprofiler_results$result)
file.name <- paste("gprof_table",i,sep = "_")
file.name <- paste(file.name,"csv",sep=".")
file.name <- paste(full.dir,file.name,sep = "/")
mydf <- apply(mydf,2,as.character)

write.csv(x=mydf,file=file.name)

if(my.terms!=FALSE){
  
  file.name <- paste("gprof",i,sep = "_")
  file.name <- paste(file.name,"pdf",sep=".")
  pdf(file = paste(full.dir,file.name,sep = "/"), width = 8,height = 6)
  # Plot for highlighting specific terms (should be run after gostplot 
  # (with interactive=FALSE))
  p <- gostplot(gprofiler_results, capped = TRUE, interactive = F)
  pp <- publish_gostplot(p, highlight_terms = my.terms, 
                         width = NA, height = NA, filename = NULL )
  pp
  dev.off()
}




