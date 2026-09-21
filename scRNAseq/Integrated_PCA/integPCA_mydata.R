library(biomaRt)
library(limma)
library(Seurat)
library(pheatmap)
library(preprocessCore)
library(corrplot)
library(dplyr)

# Run from scRNAseq/Integrated_PCA/ (or setwd accordingly)

# Load the bulk data from mESCs and NPCs
bulk.diff <- read.table(Sys.getenv("LBROMICS_BULK_DIFF_COUNTS", unset = "../../bulk/raw_counts_NPC.txt"), header = TRUE)
bulk.diff <- bulk.diff[,2:8]
bulk.diff <- as.data.frame(bulk.diff %>% 
  group_by(gene_name) %>% 
  summarise_all(funs(sum)))
rownames(bulk.diff) <- bulk.diff$gene_name
bulk.diff <- bulk.diff[,2:7]
head(bulk.diff)

bulk.und <- read.table(Sys.getenv("LBROMICS_BULK_UND_COUNTS", unset = "../../bulk/raw_counts_mESC.txt"), header = TRUE)
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

asplit <- split(1:nrow(data.combined@meta.data), data.combined@meta.data$cond)
tmp <- do.call(cbind, lapply(asplit, function(x) Matrix::rowSums(data.combined@assays$RNA@counts[,x])))
tmp <- tmp[which(rowSums(tmp) > 0),]
insilico <- tmp
cond_insilico <- as.character(data.combined@meta.data$grouping[match(colnames(insilico), data.combined@meta.data$cond)])
# Show correspondence between bulk and insilico ####
ok <- intersect(rownames(bulk.diff), rownames(insilico))
ok <- intersect(ok, rownames(bulk.und))
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

# Define sample attributes ####
batch <- c(rep("bulk mESC", ncol(bulk.und_ok)),rep("bulk NPC", ncol(bulk.diff_ok)), rep("in silico NPC", ncol(insilico_ok)))

head(insilico_ok)

cond_bulk.und <- c(rep("Lbr NT-KO", 3), rep("WT", 3))
cond_bulk.diff <- c(rep("Lbr NT-KO", 3), rep("WT", 3))
cond_insilico <- c("WT","Lbr NT-KO")

cond <- c(cond_bulk.und, cond_bulk.diff, cond_insilico)

farben <- c('#F8766D','#00BFC4')
names(farben) <- c('WT', 'Lbr NT-KO')

farben[cond]

shape <- c(1, 2, 3)
names(shape) <- c("bulk mESC", "bulk NPC", "in silico NPC")

shape[batch]

pca <- prcomp(t(merged))
pca$x[,1:2]
# Generate Fig 3c ####
pdf("./integrated_PCA.pdf",width = 8,height = 4)
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
