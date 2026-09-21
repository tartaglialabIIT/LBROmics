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
                                                         ,levels=c(0,1,2,3,4,5))
table(data.combined$integrated_snn_res.0.3)
Idents(data.combined) <- data.combined$integrated_snn_res.0.3
Idents(data.combined)

asplit <- split(1:nrow(data.combined@meta.data), data.combined@meta.data$integrated_snn_res.0.3)
tmp <- do.call(cbind, lapply(asplit, function(x) Matrix::rowSums(data.combined@assays$RNA@counts[,x])))
tmp <- tmp[which(rowSums(tmp) > 0),]
insilico <- tmp
colnames(insilico)

# Show correspondence between bulk and insilico ####
ok <- intersect(rownames(bulk.diff),rownames(bulk.und))
ok <- intersect(ok,rownames(insilico))

#ok <- intersect(ok, rownames(bulk.und))
bulk.diff_ok <- bulk.diff[ok,]
bulk.und_ok <- bulk.und[ok,]
insilico_ok <- insilico[ok, ]

#rownames(bulk_ok) <- rownames(protein_ok) <- rownames(insilico_ok) <- ok

bulk.diff_ok <- voom(bulk.diff_ok)$E
bulk.und_ok <- voom(bulk.und_ok)$E
insilico_ok <- voom(insilico_ok)$E

merged <- cbind(bulk.diff_ok, bulk.und_ok, insilico_ok)
merged <- normalize.quantiles(merged)
rownames(merged) <- rownames(bulk.diff_ok)

colnames(insilico_ok)

# Define sample attributes ####
batch <- c(rep("bulk mESC", ncol(bulk.und_ok)),
           rep("bulk NPC", ncol(bulk.diff_ok)),
           "in silico NPC_0","in silico NPC_1",
           "in silico NPC_2","in silico NPC_3",
           "in silico NPC_4","in silico NPC_5")



#farben <- c('#F8766D','#00BFC4')
#names(farben) <- c('WT', 'Lbr NT-KO')

#farben[cond]

shape <- c(1, 2, 3, 4, 5, 6, 7, 8)
names(shape) <- c("bulk mESC", "bulk NPC", "in silico NPC_0",
                  "in silico NPC_1","in silico NPC_2",
                  "in silico NPC_3","in silico NPC_4",
                  "in silico NPC_5")

shape[batch]

pca <- prcomp(t(merged))
pca$x[,1:2]
# Generate Fig 3c ####
pdf("./integrated_PCA_bycluster.pdf",width = 12,height = 8)
par(mfrow = c(1, 2))
#plot(pca$x[,1:2], col = farben[cond])
plot(pca$x[,1:2],  pch = shape[batch])
legend("topright", c(names(shape)), pch = c(shape), bty = "n", col = "black")
abline(h = 0)
abline(v = 0)
plot(pca$x[,2:3],  pch = shape[batch])
legend("topright", c(names(shape)), pch = c(shape), bty = "n", col = "black")
abline(h = 0)
abline(v = 0)
dev.off()
