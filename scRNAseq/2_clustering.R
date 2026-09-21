#!/usr/bin/env Rscript
# Working directory: set LBROMICS_SCRNA_ROOT to the Seurat analysis directory
# containing files/, RDS_objects/, etc. Falls back to getwd().
wd <- Sys.getenv("LBROMICS_SCRNA_ROOT", unset = "")
if (!nzchar(wd)) wd <- getwd()
setwd(wd)

library(Seurat)
library(ggplot2)
library(cowplot)
library(mcclust)
library(clustree)

cl2_path <- Sys.getenv("LBROMICS_GOOD_CELLS", unset = "")
if (!nzchar(cl2_path)) {
  candidates <- c(
    "./files/good_cells.csv",
    "good_cells.csv",
    file.path(dirname(wd), "good_cells.csv")
  )
  # If running from repo scRNAseq/ with deposited file:
  script_candidates <- candidates
  hit <- script_candidates[file.exists(script_candidates)]
  if (length(hit) < 1) {
    stop("good_cells.csv not found; set LBROMICS_GOOD_CELLS or place file at files/good_cells.csv")
  }
  cl2_path <- hit[[1]]
}
cl2 <- read.csv(cl2_path, row.names = 1)

raw_wt <- Sys.getenv("LBROMICS_SCRNA_WT_10X", unset = "")
raw_mut <- Sys.getenv("LBROMICS_SCRNA_MUT_10X", unset = "")
if (!nzchar(raw_wt)) {
  candidates <- c(
    file.path(dirname(wd), "raw_data", "s1_DR1"),
    file.path(wd, "..", "raw_data", "s1_DR1"),
    "raw_data/s1_DR1"
  )
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) < 1) stop("WT 10X directory not found; set LBROMICS_SCRNA_WT_10X")
  raw_wt <- hit[[1]]
}
if (!nzchar(raw_mut)) {
  candidates <- c(
    file.path(dirname(wd), "raw_data", "s2_A8"),
    file.path(wd, "..", "raw_data", "s2_A8"),
    "raw_data/s2_A8"
  )
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) < 1) stop("Mutant 10X directory not found; set LBROMICS_SCRNA_MUT_10X")
  raw_mut <- hit[[1]]
}

expression_matrix.WT <- Read10X(data.dir = raw_wt)
seurat_object.WT <- CreateSeuratObject(counts = expression_matrix.WT)
seurat_object.WT$cond <- "WT"

colnames(seurat_object.WT)
# Select only good cells
seurat_object.WT[["CellName"]] <- paste(colnames(seurat_object.WT), "WT", sep="-") 
seurat_object.WT <- subset(seurat_object.WT, subset = CellName %in% rownames(cl2) )

seurat_object.WT <- subset(seurat_object.WT, subset = nFeature_RNA > 0)
seurat_object.WT <- NormalizeData(seurat_object.WT, verbose = FALSE)
seurat_object.WT <- FindVariableFeatures(seurat_object.WT, selection.method = "vst", nfeatures = 2000)

expression_matrix.mutant <- Read10X(data.dir = raw_mut)
seurat_object.mutant <- CreateSeuratObject(counts = expression_matrix.mutant)
seurat_object.mutant$cond <- "mutant"
# Select only good cells
seurat_object.mutant[["CellName"]] <- paste(colnames(seurat_object.mutant), "mutant", sep="-") 
seurat_object.mutant <- subset(seurat_object.mutant, subset = CellName %in% rownames(cl2) )

# Select only good cells
seurat_object.mutant <- subset(seurat_object.mutant, subset = nFeature_RNA > 0)
seurat_object.mutant <- NormalizeData(seurat_object.mutant, verbose = FALSE)
seurat_object.mutant <- FindVariableFeatures(seurat_object.mutant, selection.method = "vst", nfeatures = 2000)

my.anchors <- FindIntegrationAnchors(object.list = list(seurat_object.WT, seurat_object.mutant),
                                          dims = 1:20, anchor.features = 2000)
data.combined <- IntegrateData(anchorset = my.anchors, dims = 1:20)
DefaultAssay(data.combined) <- "integrated"
data.combined <- ScaleData(data.combined, verbose = FALSE)

data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:20,set.seed=42)
data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:20)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:20,k.param=20)

# Find clusters at different resolutions
data.combined <- FindClusters(data.combined, resolution =0.3)

pdf(file = "./plots/clustree.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
clustree(data.combined,prefix = 'integrated_snn_res.' )
dev.off()

# Use clustree to visualize the clusterings
table(Idents(data.combined))
table(Idents(data.combined), data.combined$cond)

data.combined$cond <- factor(data.combined$cond, levels = c('WT', 'mutant'))
data.combined@meta.data$integrated_snn_res.0.2 <- factor(x = data.combined@meta.data$integrated_snn_res.0.2
                                                         ,levels=c(0,2,3,1))
data.combined@meta.data$integrated_snn_res.0.3 <- factor(x = data.combined@meta.data$integrated_snn_res.0.3
                                                         ,levels=c(1,2,3,6,4,0,5))

library(repr)
data.combined[["percent.mt"]] <- PercentageFeatureSet(data.combined, pattern = "^mt-", assay = 'RNA')
data.combined[["percent.ribo"]] <- PercentageFeatureSet(data.combined, pattern = "^Rp[sl]", assay = 'RNA')

data.combined$nFeature_RNA

pdf(file = "./plots_JULY_2023/violin_percent_mt_percent_ribo.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
p1 <- VlnPlot(data.combined,features = c("percent.mt","percent.ribo"),
              group.by = "integrated_snn_res.0.3",pt.size = 0.1,split.by = "cond")
p1
dev.off()
pdf(file = "./plots_JULY_2023/violin_nCountRNA_nFeature_RNA.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
p1 <- VlnPlot(data.combined,features = c("nCount_RNA","nFeature_RNA"),
              group.by = "integrated_snn_res.0.3",pt.size = 0.1,split.by = "cond")
p1
dev.off()

# We remove cells belonging to cluster 6
Idents(data.combined)
data.combined <- SetIdent(data.combined, value = "integrated_snn_res.0.3")
Idents(data.combined)
sub_data.combined <- subset(x = data.combined, idents = 6, invert = TRUE)
saveRDS(sub_data.combined,"./RDS_objects/sub_data_combined_and_clustered.rds")

Idents(sub_data.combined)

pdf(file = "./plots_JULY_2023/umap_res01.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1 <- DimPlot(sub_data.combined, reduction = "umap", group.by = "integrated_snn_res.0.1")
p1
dev.off()
pdf(file = "./plots_JULY_2023/umap_res02.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1 <- DimPlot(sub_data.combined, reduction = "umap", group.by = "integrated_snn_res.0.2")
p1
dev.off()
pdf(file = "./plots_JULY_2023/umap_res03.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1 <- DimPlot(sub_data.combined, reduction = "umap", group.by = "integrated_snn_res.0.3")
p1
dev.off()
pdf(file = "./plots_JULY_2023/umap_cond.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1 <- DimPlot(sub_data.combined, reduction = "umap", group.by = "cond")
p1
dev.off()
pdf(file = "./plots_JULY_2023/umap_xist.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1 <- FeaturePlot(sub_data.combined, features = 'Xist' )
p1
dev.off()
