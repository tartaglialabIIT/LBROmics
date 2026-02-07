### make sure you have WebGestaltR installed in your system
library("WebGestaltR")
library("dplyr")

options(bitmapType='cairo')

### creates output directories if they do not exist
make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
} 

### I wrapped webgestaltR into testfunction to handle possible errors
testFunction <- function (...) {
  return(tryCatch(WebGestaltR(...), error=function(e) NULL))
}

### this function launches webgestalt using the interesting gene list (iList), the background (universe). projName is the name of the output directory, isOut is set to true if
### results have to be written in files. If gene IDs are ensembl gene ids with version number, this number is removed
### functional enrichment analysis will be performed on all the categories in CATEGORIES global variable, looping over such categories and the corresponding output files
### stored in OUT.DIR.ORA

goEnrich <- function(iList,universe,projName,isOut,...){
  if (isOut==TRUE) {
    lapply(paste(OUT.DIR.ORA,"Project_",projName,sep=""),unlink,recursive=TRUE )
  }
  iList = sub(iList,pattern = "\\.\\d+$", replacement = "")
  universe = sub(universe,pattern = "\\.\\d+$", replacement = "")

  wg = mapply(testFunction, enrichDatabase = CATEGORIES , outputDirectory = OUT.DIR.ORA,
              MoreArgs=list(enrichMethod="ORA", organism="mmusculus",interestGene=iList, # sigMethod = "top",
                            interestGeneType="ensembl_gene_id",referenceGene=universe,referenceGeneType="ensembl_gene_id", 
                            minNum = 5, maxNum=2000, nThreads=20, isOutput=isOut, projectName=projName,...) )
  return(wg)
}

### this function, similar to goEnrich, performs GSEA using the log2FoldChange as ranking metric (make sure this has been shrunk with apglm) 

gseaEnrich <- function(annotDE,projName,...){
  lapply(paste(OUT.DIR.GSEA,"Project_",projName,sep=""),unlink,recursive=TRUE )
  rankList = data.frame(Genes=row.names(annotDE),
                        #DDE = -log10(annotDE$pvalue) * sign(annotDE$log2FoldChange),
                        DDE = annotDE$log2FoldChange,
                        stringsAsFactors = FALSE)
  
  rankList$Genes = sub(rankList$Genes,pattern = "\\.\\d+$", replacement = "")
  wg = mapply(testFunction, enrichDatabase = CATEGORIES , outputDirectory = OUT.DIR.GSEA,
              MoreArgs=list(enrichMethod="GSEA", organism="mmusculus",interestGene=rankList, 
                            interestGeneType="ensembl_gene_id",  minNum = 10, maxNum=1000, reportNum = 40, isOutput=TRUE,
                            nThreads=20, projectName=projName,... ))
  return(wg)
}

### launch ORA analysis on all the differentially expressed genes, only the UPregulated and only the DOWNregulated. resAnno is a deseq2 output table which also
### contains a DEREG_FLAG folder specifying if the gene is upregulated ("UP"), downregulated ("DOWN") or not affected ("NO")
### projNAME is a string describing the pairwise comparison for which DE analysis was performed (e.g. "WT_vs_KO")

oraLaunch <- function(resAnno,projNAME, fdrTHR) {
  oraRes <- goEnrich(rownames(subset(resAnno, DEREG_FLAG != "NO"  )),
                        rownames(resAnno),
                     projNAME,
                        TRUE,
                        fdrThr = fdrTHR
  )
  
  oraResUp <- goEnrich(rownames(subset(resAnno, DEREG_FLAG == "UP"  )),
                           rownames(resAnno),
                           paste(projNAME,"UP",sep="_"),
                           TRUE,
                           fdrThr = fdrTHR
  )
  
  oraResDown <- goEnrich(rownames(subset(resAnno, DEREG_FLAG == "DOWN"  )),
                             rownames(resAnno),
                             paste(projNAME,"DOWN",sep="_"),
                             TRUE,
                             fdrThr = fdrTHR
  )
  return(list(oraRes,oraResUp,oraResDown))
}

### mouse functional categories
#CATEGORIES = c("geneontology_Biological_Process",
#               "chromosomalLocation_CytogeneticBand")

CATEGORIES = c("geneontology_Biological_Process","geneontology_Biological_Process_noRedundant",
              "geneontology_Cellular_Component",
              "geneontology_Cellular_Component_noRedundant",
              "geneontology_Molecular_Function",
              "geneontology_Molecular_Function_noRedundant",
              "pathway_KEGG",
              "pathway_Panther",
              "pathway_Reactome",
              "network_CORUM",
              "pathway_Wikipathway",
              "network_PPI_BIOGRID",
              "network_Transcription_Factor_target",
              "chromosomalLocation_CytogeneticBand")
### before launching the script, create a deseq2 folder (usually i put in ../files/deseq2/) in which I put the "GO/webgestalt/" folder, which in turns contain the
### ORA and GSEA folders

dir.deseq = "/mnt/large/jfiorentino/Cerase_Data/deseq2"
dir.create(dir.deseq)


### specify the ORA folder and create one folder for each category. When you will perform ORA analysis, the result folder/files for each category will be put
### in the corresponding folders. For each project (differential expression comparison) you will have a different folder under the same category folder

DIR.ORA <- paste(dir.deseq,"/GO/webgestalt/ORA/",sep="")
OUT.DIR.ORA <- mapply(paste, DIR.ORA,CATEGORIES,"/", MoreArgs=list(sep="") , SIMPLIFY=TRUE)
mapply(make.dir, OUT.DIR.ORA) 

### specify the ORA folder and create one folder for each category

DIR.GSEA <- paste(dir.deseq,"/GO/webgestalt/GSEA/",sep="")
OUT.DIR.GSEA <- mapply(paste, DIR.GSEA,CATEGORIES,"/", MoreArgs=list(sep="") , SIMPLIFY=TRUE)
mapply(make.dir, OUT.DIR.GSEA) 

### here you put the code for the differential expression. At the end, you have to produce, for each comparison, a table containing the DE results with column that
### annotates the genes based on the differential expression status (see oraLaunch documentation). 
#wt_mutant.annotated.sc <-read.csv("/mnt/large/jfiorentino/Cerase_Data/scRNA_seq/in_silico_bulk_DE/res_tot_DESeq2_scparams_with_IDs.csv")
#print(head(wt_mutant.annotated.sc))
#wt_mutant.annotated.sc$DEREG_FLAG <- "NO"
#wt_mutant.annotated.sc[which(wt_mutant.annotated.sc$log2FoldChange>0.58 & wt_mutant.annotated.sc$padj<0.05),]$DEREG_FLAG <- "UP"
#wt_mutant.annotated.sc[which(wt_mutant.annotated.sc$log2FoldChange< -0.58 & wt_mutant.annotated.sc$padj<0.05),]$DEREG_FLAG  <- "DOWN"

#print(colnames(wt_mutant.annotated.sc))
#print(head(wt_mutant.annotated.sc))
#print(wt_mutant.annotated.sc[wt_mutant.annotated.sc[,"Ensembl.gene.ID"]=='NaN',])

#wt_mutant.annotated.sc <- wt_mutant.annotated.sc["ENSMUS" %in% wt_mutant.annotated.sc[,"Ensembl.gene.ID"],]
#wt_mutant.annotated.sc <- wt_mutant.annotated.sc[!(is.na(wt_mutant.annotated.sc$Ensembl.gene.ID) | wt_mutant.annotated.sc$Ensembl.gene.ID==""), ]

#row.names(wt_mutant.annotated.sc) <- wt_mutant.annotated.sc[,'Ensembl.gene.ID']
#print(length(rownames(wt_mutant.annotated.sc)))
### Now we perform ORA analysis. Suppose such table is stored in the wt_ko.annotated object,
### with project name "WT_vs_KO" (you will have to repeat this for each comparison). oraLaunch will perform ORA analysis and store the result tables in a list
### It will also store the result files in the specified directories
### the first element contains the result table for all the differentially expressed genes, and will be assigned to ora.all variable.
### the second element contains the result table for the upregulated genes, and will be assigned to ora.up variable.
### the third element contains the result table for the downregulated genes, and will be assigned to ora.down variable.
### if you have multiple deseq2 tables, make sure that these variables have with different names

#ora.results.sc = oraLaunch(wt_mutant.annotated.sc,"WT_vs_MUTANT_sc",0.05) ### 0.05 is the FDR threshold for the enriched categories. Set it as you like, I recommend a higher values, since you can always filter afterwards
#ora.all.sc <- ora.results.sc[[1]]
#ora.up.sc <- ora.results.sc[[2]]
#ora.down.sc <- ora.results.sc[[3]]

### Now we do GSEA analysis

### now we perform the GSEA analysis. We will use the shrunk log2foldchange as a ranking metric. instead of filtering based on fdr, we report the top 40 categories.
### Again, these are reported in the gsea.results tables, as well as in the output folders relative to each category

#gsea.results.sc <- gseaEnrich(wt_mutant.annotated.sc,"WT_vs_MUTANT_sc",sigMethod = "top",  topThr = 40)




### here you put the code for the differential expression. At the end, you have to produce, for each comparison, a table containing the DE results with column that
### annotates the genes based on the differential expression status (see oraLaunch documentation). 
wt_mutant.annotated.und <-read.csv("/mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/DE_analysis_results/deseq2_results_mESC.txt",sep = '\t')
wt_mutant.annotated.und$DEREG_FLAG <- "NO"
wt_mutant.annotated.und[which(wt_mutant.annotated.und$log2FoldChange>1.0 & wt_mutant.annotated.und$padj<0.01),]$DEREG_FLAG <- "UP"
wt_mutant.annotated.und[which(wt_mutant.annotated.und$log2FoldChange< -1.0 & wt_mutant.annotated.und$padj<0.01),]$DEREG_FLAG  <- "DOWN"
row.names(wt_mutant.annotated.und) <- wt_mutant.annotated.und$gene_id
### Now we perform ORA analysis. Suppose such table is stored in the wt_ko.annotated object,
### with project name "WT_vs_KO" (you will have to repeat this for each comparison). oraLaunch will perform ORA analysis and store the result tables in a list
### It will also store the result files in the specified directories
### the first element contains the result table for all the differentially expressed genes, and will be assigned to ora.all variable.
### the second element contains the result table for the upregulated genes, and will be assigned to ora.up variable.
### the third element contains the result table for the downregulated genes, and will be assigned to ora.down variable.
### if you have multiple deseq2 tables, make sure that these variables have with different names

ora.results.und = oraLaunch(wt_mutant.annotated.und,"WT_vs_MUTANT_bulk_und",0.05) ### 0.05 is the FDR threshold for the enriched categories. Set it as you like, I recommend a higher values, since you can always filter afterwards
ora.all.und <- ora.results.und[[1]]
ora.up.und <- ora.results.und[[2]]
ora.down.und <- ora.results.und[[3]]

### Now we do GSEA analysis

### now we perform the GSEA analysis. We will use the shrunk log2foldchange as a ranking metric. instead of filtering based on fdr, we report the top 40 categories.
### Again, these are reported in the gsea.results tables, as well as in the output folders relative to each category

gsea.results.und <- gseaEnrich(wt_mutant.annotated.und,"WT_vs_MUTANT_bulk_und",sigMethod = "top",  topThr = 40)


### here you put the code for the differential expression. At the end, you have to produce, for each comparison, a table containing the DE results with column that
### annotates the genes based on the differential expression status (see oraLaunch documentation). 
wt_mutant.annotated.diff <- read.csv("/mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/DE_analysis_results/deseq2_results_NPC.txt",sep = '\t')
wt_mutant.annotated.diff$DEREG_FLAG <- "NO"
wt_mutant.annotated.diff[which(wt_mutant.annotated.diff$log2FoldChange>1.0 & wt_mutant.annotated.diff$padj<0.01),]$DEREG_FLAG <- "UP"
wt_mutant.annotated.diff[which(wt_mutant.annotated.diff$log2FoldChange< -1.0 & wt_mutant.annotated.diff$padj<0.01),]$DEREG_FLAG  <- "DOWN"
row.names(wt_mutant.annotated.diff) <- wt_mutant.annotated.diff$gene_id
### Now we perform ORA analysis. Suppose such table is stored in the wt_ko.annotated object,
### with project name "WT_vs_KO" (you will have to repeat this for each comparison). oraLaunch will perform ORA analysis and store the result tables in a list
### It will also store the result files in the specified directories
### the first element contains the result table for all the differentially expressed genes, and will be assigned to ora.all variable.
### the second element contains the result table for the upregulated genes, and will be assigned to ora.up variable.
### the third element contains the result table for the downregulated genes, and will be assigned to ora.down variable.
### if you have multiple deseq2 tables, make sure that these variables have with different names

ora.results.diff = oraLaunch(wt_mutant.annotated.diff,"WT_vs_MUTANT_bulk_diff",0.05) ### 0.05 is the FDR threshold for the enriched categories. Set it as you like, I recommend a higher values, since you can always filter afterwards
ora.all.diff <- ora.results.diff[[1]]
ora.up.diff <- ora.results.diff[[2]]
ora.down.diff <- ora.results.diff[[3]]

### Now we do GSEA analysis

### now we perform the GSEA analysis. We will use the shrunk log2foldchange as a ranking metric. instead of filtering based on fdr, we report the top 40 categories.
### Again, these are reported in the gsea.results tables, as well as in the output folders relative to each category

gsea.results.diff <- gseaEnrich(wt_mutant.annotated.diff,"WT_vs_MUTANT_bulk_diff",sigMethod = "top",  topThr = 40)
