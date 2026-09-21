library(biomaRt)
library(limma)
library(Seurat)
library(pheatmap)
library(preprocessCore)
library(corrplot)
library(dplyr)

InsilicoCl <- function(data,cl){
  data.tmp <- subset(x=data,idents = cl)
  asplit <- split(1:nrow(data.tmp@meta.data), data.tmp@meta.data$cond)
  tmp <- do.call(cbind, lapply(asplit, function(x) Matrix::rowSums(data.tmp@assays$RNA@counts[,x])))
  tmp <- tmp[which(rowSums(tmp) > 0),]
  insilico <- tmp
  colnames(insilico)
  
  
  #cond_insilico <- as.character(data.tmp@meta.data$grouping[match(colnames(insilico), data.tmp@meta.data$cond)])
  insilico
}

wd <- "/Users/jonathan/Desktop/IIT/Cerase_single_cell/ANALYSIS/Seurat_analysis/Integrated_PCA/"
setwd(wd)

# Load the bulk data from mESCs and NPCs
bulk.diff <- read.table("../../../bulk_RNA_seq/analysis_results/raw_counts/raw_counts_diff.txt", header = T)
bulk.diff <- bulk.diff[,2:8]
bulk.diff <- as.data.frame(bulk.diff %>% 
  group_by(gene_name) %>% 
  summarise_all(funs(sum)))
rownames(bulk.diff) <- bulk.diff$gene_name
bulk.diff <- bulk.diff[,2:7]
head(bulk.diff)

bulk.und <- read.table("../../../bulk_RNA_seq/analysis_results/raw_counts/raw_counts_und.txt", header = T)
bulk.und <- bulk.und[,2:8]
bulk.und <- as.data.frame(bulk.und %>% 
                             group_by(gene_name) %>% 
                             summarise_all(funs(sum)))
rownames(bulk.und) <- bulk.und$gene_name
bulk.und <- bulk.und[,2:7]
head(bulk.und)

cl2 <- read.csv("../files/leiden_clusters.csv",row.names = 1)


# Load and generate in silico bulk from scRNA-seq data
data.combined <- readRDS("../RDS_objects/sub_data_combined_and_clustered.rds")
data.combined@meta.data$integrated_snn_res.0.3 <- factor(x = data.combined@meta.data$integrated_snn_res.0.3
                                                         ,levels=c(1,2,3,4,0,5))
table(data.combined$integrated_snn_res.0.3)
Idents(data.combined) <- data.combined$integrated_snn_res.0.3
Idents(data.combined)

insilico0 <- InsilicoCl(data.combined,0)
insilico1 <- InsilicoCl(data.combined,1)
insilico2 <- InsilicoCl(data.combined,2)
insilico3 <- InsilicoCl(data.combined,3)
insilico4 <- InsilicoCl(data.combined,4)
insilico5 <- InsilicoCl(data.combined,5)

# Show correspondence between bulk and insilico ####
ok <- intersect(rownames(bulk.diff),rownames(bulk.und))
ok <- intersect(ok,rownames(insilico0))
ok <- intersect(ok,rownames(insilico1))
ok <- intersect(ok,rownames(insilico2))
ok <- intersect(ok,rownames(insilico3))
ok <- intersect(ok,rownames(insilico4))
ok <- intersect(ok,rownames(insilico5))

#ok <- intersect(ok, rownames(bulk.und))
bulk.diff_ok <- bulk.diff[ok,]
bulk.und_ok <- bulk.und[ok,]
insilico_ok0 <- insilico0[ok, ]
insilico_ok1 <- insilico1[ok, ]
insilico_ok2 <- insilico2[ok, ]
insilico_ok3 <- insilico3[ok, ]
insilico_ok4 <- insilico4[ok, ]
insilico_ok5 <- insilico5[ok, ]

#rownames(bulk_ok) <- rownames(protein_ok) <- rownames(insilico_ok) <- ok

bulk.diff_ok <- voom(bulk.diff_ok)$E
bulk.und_ok <- voom(bulk.und_ok)$E
insilico_ok0 <- voom(insilico_ok0)$E
insilico_ok1 <- voom(insilico_ok1)$E
insilico_ok2 <- voom(insilico_ok2)$E
insilico_ok3 <- voom(insilico_ok3)$E
insilico_ok4 <- voom(insilico_ok4)$E
insilico_ok5 <- voom(insilico_ok5)$E
merged <- cbind(bulk.diff_ok, bulk.und_ok, insilico_ok0,
                insilico_ok1, insilico_ok2, insilico_ok3,
                insilico_ok4, insilico_ok5)
merged <- normalize.quantiles(merged)
rownames(merged) <- rownames(bulk.diff_ok)

# Define sample attributes ####
batch <- c(rep("bulk mESC", ncol(bulk.und_ok)),
           rep("bulk NPC", ncol(bulk.diff_ok)), 
           rep("in silico NPC_0", ncol(insilico_ok0)), 
           rep("in silico NPC_1", ncol(insilico_ok1)), 
           rep("in silico NPC_2", ncol(insilico_ok2)), 
           rep("in silico NPC_3", ncol(insilico_ok3)), 
           rep("in silico NPC_4", ncol(insilico_ok4)), 
           rep("in silico NPC_5", ncol(insilico_ok5)))

head(insilico_ok0)
head(bulk.und_ok)
head(bulk.diff_ok)

table(data.combined$cond)

cond_bulk.und <- c(rep("Lbr NT-KO", 3), rep("WT", 3))
cond_bulk.diff <- c(rep("Lbr NT-KO", 3), rep("WT", 3))
cond_insilico0 <- c("WT","Lbr NT-KO")
cond_insilico1 <- c("WT","Lbr NT-KO")
cond_insilico2 <- c("WT","Lbr NT-KO")
cond_insilico3 <- c("WT","Lbr NT-KO")
cond_insilico4 <- c("WT","Lbr NT-KO")
cond_insilico5 <- c("WT","Lbr NT-KO")

cond <- c(cond_bulk.und, cond_bulk.diff, cond_insilico0,
          cond_insilico1,cond_insilico2,cond_insilico3,
          cond_insilico4,cond_insilico5)

farben <- c('#F8766D','#00BFC4')
names(farben) <- c('WT', 'Lbr NT-KO')

farben[cond]

shape <- c(1, 2, 3, 4, 5, 6, 7, 8)
names(shape) <- c("bulk mESC", "bulk NPC", "in silico NPC_0",
                  "in silico NPC_1","in silico NPC_2",
                  "in silico NPC_3","in silico NPC_4",
                  "in silico NPC_5")

shape[batch]

pca <- prcomp(t(merged))
pca$x[,1:2]
# Generate Fig 3c ####
pdf("./integrated_PCA_bycluster_and_cond.pdf",width = 12,height = 8)
par(mfrow = c(1, 2))
#plot(pca$x[,1:2], col = farben[cond])
plot(pca$x[,1:2], col = farben[cond], pch = shape[batch])
legend("topright", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c('#F8766D','#00BFC4', "black", "black", "black"))
abline(h = 0)
abline(v = 0)
plot(pca$x[,2:3], col = farben[cond], pch = shape[batch])
legend("topright", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c('#F8766D','#00BFC4', "black", "black", "black"))
abline(h = 0)
abline(v = 0)
dev.off()
