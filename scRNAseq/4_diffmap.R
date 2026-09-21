#!/usr/bin/env Rscript
library(Seurat) 
library(SingleCellExperiment)
library(destiny)

# Set LBROMICS_SCRNA_ROOT to the Seurat analysis directory.
wd <- Sys.getenv("LBROMICS_SCRNA_ROOT", unset = "")
if (!nzchar(wd)) wd <- getwd()
# Optional alternate working directory historically used for destiny runs:
# wd <- Sys.getenv("LBROMICS_DIFFMAP_DIR", unset = wd)

setwd(wd)

seuratObject <- readRDS("./RDS_objects/sub_data_combined_and_clustered.rds")
seuratObject@meta.data$integrated_snn_res.0.3 <- factor(x = seuratObject@meta.data$integrated_snn_res.0.3
                                                         ,levels=c(1,2,3,4,0,5))
Idents(seuratObject) <- seuratObject@meta.data$integrated_snn_res.0.3

sce <- as.SingleCellExperiment(seuratObject)
print(sce)
#this has the cell classification
print(table(sce$ident))
dm <- DiffusionMap(sce, verbose = TRUE)

save(dm, file = "./files/DiffusionMap.RData")

tmp <- data.matrix(data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2],DC3 = eigenvectors(dm)[, 3], row.names = colnames(seuratObject)))
seuratObject[["dm"]] <- CreateDimReducObject(embeddings = tmp, key="DC_", assay=DefaultAssay(seuratObject))
saveRDS(seuratObject,"./files/sub_data_combined_and_clustered_dm.rds")

load("./files/DiffusionMap.RData")
dm
library(ggplot2)
cellLabels <- sce$ident
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Samples = cellLabels)
pdf("./plots_JULY_2023/trajectory_inference/DC1_DC2.pdf", w=11, h=8.5)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Samples)) +
  geom_point()  + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()
dev.off()