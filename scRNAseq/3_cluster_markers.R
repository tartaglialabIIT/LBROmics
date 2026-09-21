#!/usr/bin/env Rscript
# Set LBROMICS_SCRNA_ROOT to the Seurat analysis directory (contains RDS_objects/).
wd <- Sys.getenv("LBROMICS_SCRNA_ROOT", unset = "")
if (!nzchar(wd)) wd <- getwd()
setwd(wd)

library(Seurat)
library(ggplot2)
library(cowplot)
library(mcclust)
library(clustree)
library(dplyr)
source("myfunctions.R")

# Load the clustered rds object
data.combined <- readRDS("./RDS_objects/sub_data_combined_and_clustered.rds")
data.combined


table(data.combined$integrated_snn_res.0.3)

data.combined$cond <- factor(data.combined$cond, levels = c('WT', 'mutant'))
data.combined@meta.data$integrated_snn_res.0.2 <- factor(x = data.combined@meta.data$integrated_snn_res.0.2
                                                         ,levels=c(0,2,3,1))
data.combined@meta.data$integrated_snn_res.0.3 <- factor(x = data.combined@meta.data$integrated_snn_res.0.3
                                                         ,levels=c(1,2,3,4,0,5))

df <- as.data.frame(table(data.combined$integrated_snn_res.0.3,data.combined$cond))
colnames(df) <- c("cluster","cond","count")
write.csv(x=df,file = "./files_JULY_2023/cond_df.csv")

sum.WT <- sum(df[df$cond=="WT","count"])
sum.mutant <- sum(df[df$cond=="mutant","count"])

df[df$cond=="WT","count"] <- df[df$cond=="WT","count"]/sum.WT
df[df$cond=="mutant","count"] <- df[df$cond=="mutant","count"]/sum.WT

df.new <- transform(df, z = count / ave(count, cluster, FUN = sum))
library(ggplot2)


library(patchwork)

# Define custom colors
custom_colors <- c(WT = "black", mutant = "#E72226")
DefaultAssay(data.combined) <- "RNA"

# Create each plot with custom colors
p1 <- VlnPlot(
  data.combined,
  features = "Dppa5a",
  group.by = "integrated_snn_res.0.3",
  pt.size = 0.0,
  split.by = "cond"
) + scale_fill_manual(values = custom_colors) + theme(legend.position = "none")

p2 <- VlnPlot(
  data.combined,
  features = "Dnmt3l",
  group.by = "integrated_snn_res.0.3",
  pt.size = 0.0,
  split.by = "cond"
) + scale_fill_manual(values = custom_colors) + theme(legend.position = "none")

p3 <- VlnPlot(
  data.combined,
  features = "Dppa3",
  group.by = "integrated_snn_res.0.3",
  pt.size = 0.0,
  split.by = "cond"
) + scale_fill_manual(values = custom_colors) + theme(legend.position = "none")

p4 <- VlnPlot(
  data.combined,
  features = "Klf5",
  group.by = "integrated_snn_res.0.3",
  pt.size = 0.0,
  split.by = "cond"
) + scale_fill_manual(values = custom_colors) + theme(legend.position = "none")


# Combine plots side by side
combined_plot <- p1 | p2 | p3 | p4

# Save to PDF
pdf(file = "./plots_JULY_2023/markers_figs6.pdf",
    width = 13, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches
print(combined_plot)
dev.off()


# Stacked
pdf("./plots_JULY_2023/norm_barplot_condition.pdf",width = 6,height = 4)
ggplot(df.new, aes(fill=cond, y=z, x=cluster)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic()+
  xlab("Cluster") +
  ylab("Normalized frequency")
dev.off()

DefaultAssay(data.combined) <- "RNA"
VlnPlot(data.combined,features = "Tpm3",split.by = "cond")
VlnPlot(data.combined,features = "AY036118",split.by = "cond")
FeaturePlot(data.combined,features = "Lbr")
DefaultAssay(data.combined) <- "RNA"


naive_ground <- c("Klf17", "Dppa5a", "Dnmt3l",
                  "Gata6", "Tbx3", "Il6st", "Dppa3", 
                  "Klf5")
pdf(file = "./plots_JULY_2023/cell_anno/final_figs/violin_naive_pluripotecy_genes.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 9) # The height of the plot in inches
p <- VlnPlot(data.combined,features = naive_ground,group.by = "integrated_snn_res.0.3",split.by = "cond")
p 
dev.off()

primed <- c("Cd24a", "Zic2", "Sfrp2","Otx2")
pdf(file = "./plots_JULY_2023/cell_anno/final_figs/violin_primed_pluripotecy_genes.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p <- VlnPlot(data.combined,features = primed,group.by = "integrated_snn_res.0.3",split.by = "cond")
p
dev.off()

neuroepi <- c("Nes","Sox2","Notch1","Hes1","Hes3","Ocln","Sox10","Cdh1")
pdf(file = "./plots_JULY_2023/cell_anno/final_figs/violin_neuroepi_genes.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p <- VlnPlot(data.combined,features = neuroepi,group.by = "integrated_snn_res.0.3",split.by = "cond")
p
dev.off()


meta.data <- data.combined[[]]

counts <- group_by(meta.data, cond, integrated_snn_res.0.1) %>% summarise(count = n())
ggplot(counts, aes(integrated_snn_res.0.1, count, fill = cond)) +
  geom_bar(stat = 'identity')
ggsave(filename = "./plots_JULY_2023/cluster_comp_res01.pdf")
counts <- group_by(meta.data, cond, integrated_snn_res.0.2) %>% summarise(count = n())
ggplot(counts, aes(integrated_snn_res.0.2, count, fill = cond)) +
  geom_bar(stat = 'identity')
ggsave(filename = "./plots_JULY_2023/cluster_comp_res02.pdf")
counts <- group_by(meta.data, cond, integrated_snn_res.0.3) %>% summarise(count = n())
ggplot(counts, aes(integrated_snn_res.0.3, count, fill = cond)) +
  geom_bar(stat = 'identity')
ggsave(filename = "./plots_JULY_2023/cluster_comp_res03.pdf")

# Compute cluster specific markers at different resolutions
markers01 <- list()
clusters <- table(data.combined[['integrated_snn_res.0.1']])
for (i in 0:(length(clusters)-1)){
  df <- FindSpecificMarkers(data.combined,"integrated_snn_res.0.1",i)
  markers01[[i+1]] <- df
}

markers02 <- list()
clusters <- table(data.combined[['integrated_snn_res.0.2']])
for (i in 0:(length(clusters)-1)){
  df <- FindSpecificMarkers(data.combined,"integrated_snn_res.0.2",i)
  markers02[[i+1]] <- df
}

table(data.combined[['integrated_snn_res.0.3']])

markers03 <- list()
clusters <- table(data.combined[['integrated_snn_res.0.3']])
for (i in 0:(length(clusters)-1)){
  df <- FindSpecificMarkers(data.combined,"integrated_snn_res.0.3",i)
  markers03[[i+1]] <- df
}



markers03

for (i in 1: length(markers01)){
  print(length(rownames(markers01[[i]])))
}
for (i in 1: length(markers02)){
  print(length(rownames(markers02[[i]])))
}
for (i in 1: length(markers03)){
  print(length(rownames(markers03[[i]])))
  print(head(markers03[[i]]))
}

library(openxlsx)
write.xlsx(markers03, "Specific_cluster_markers.xlsx")

DefaultAssay(data.combined) <- "RNA"
data.combined <- ScaleData(data.combined)
# Make a heatmap of the Xlr genes
row.names(data.combined)[grep("Xlr", row.names(data.combined))]

"Xlr" %in% row.names(data.combined)
pdf(file = "./plots_JULY_2023/dotplot_Xlr.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches
p1 <- DotPlot(object = data.combined, features = row.names(data.combined)[grep("Xlr", row.names(data.combined))],
        group.by = "cond")
p1
dev.off()
DoHeatmap(object = data.combined,features =row.names(data.combined)[grep("Xlr", row.names(data.combined))],
          group.by = "cond")
# Make a heatmap of NPC markers 
NPC_markers <- c('Pax3','Pax6','Sox1','Sox2','Otx2','Ascl1','Smarca4','Msi1','Msi2','Nes')
pdf(file = "./plots_JULY_2023/dotplot_NPC_markers_cond.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches
p1 <- DotPlot(object = data.combined, features = intersect(row.names(data.combined),NPC_markers),
              group.by = "cond")
p1
dev.off()
pdf(file = "./plots_JULY_2023/dotplot_NPC_markers.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches
p1 <- DotPlot(object = data.combined, features = intersect(row.names(data.combined),NPC_markers),
              group.by = "integrated_snn_res.0.3")
p1
dev.off()


# Make heatmaps of cluster markers at different resolutions
data.combined@meta.data$integrated_snn_res.0.1 <- factor(x = data.combined@meta.data$integrated_snn_res.0.1
                                                     ,levels=c(0,1,2))


data.combined$cond <- factor(data.combined$cond, levels = c('WT', 'mutant'))



my.features01 <- c(markers01[[1]]$Gene[1:10],markers01[[2]]$Gene[1:10],markers01[[3]]$Gene[1:10])
DoMultiBarHeatmap(data.combined, features=my.features01, group.by='integrated_snn_res.0.1',label=F,
                  assay='RNA',slot='scale.data',disp.max=2.5,
                  additional.group.by = c('cond'), additional.group.sort.by = c('cond'))
ggsave(filename="./plots/specmarkheatmap_res01.pdf", width = 8, height = 8, units = "in")

length(markers03)

data.combined@meta.data$integrated_snn_res.0.2 <- factor(x = data.combined@meta.data$integrated_snn_res.0.2
                                                         ,levels=c(0,2,3,1))
my.features02 <- c(markers02[[1]]$Gene[1:10],markers02[[3]]$Gene[1:10],markers02[[4]]$Gene[1:10],markers02[[2]]$Gene[1:10])
DoMultiBarHeatmap(data.combined, features=my.features02, group.by='integrated_snn_res.0.2',label=F,
                  assay='RNA',slot='scale.data',disp.max=2.5,
                  additional.group.by = c('cond'), additional.group.sort.by = c('cond'))
ggsave(filename="./plots/specmarkheatmap_res02.pdf", width = 8, height = 11, units = "in")

data.combined@meta.data$integrated_snn_res.0.3 <- factor(x = data.combined@meta.data$integrated_snn_res.0.3
                                                         ,levels=c(1,2,3,4,0,5))
my.features03 <- c(markers03[[2]]$Gene[1:5],markers03[[3]]$Gene[1:6],markers03[[4]]$Gene[1:6],
                   markers03[[5]]$Gene[1:5],markers03[[1]]$Gene[1:5],markers03[[6]]$Gene[1:5])
DoMultiBarHeatmap(data.combined, features=my.features03, group.by='integrated_snn_res.0.3',label=F,
                  assay='RNA',slot='scale.data',disp.max=2.5,
                  additional.group.by = c('cond'), additional.group.sort.by = c('cond')) + 
  theme(text = element_text(size = 14))
ggsave(filename="./plots_JULY_2023/specmarkheatmap_res03.pdf", width = 8, height = 9, units = "in")


# GO of the specific marker genes for each cluster (and each clustering resolution)
clusters <- table(data.combined[['integrated_snn_res.0.1']])
for (i in 0:(length(clusters)-1)){
  custom.GO(markers01[[i+1]]$Gene,i,"res01",rownames(data.combined))
}
clusters <- table(data.combined[['integrated_snn_res.0.2']])
for (i in 0:(length(clusters)-1)){
  custom.GO(markers02[[i+1]]$Gene,i,"res02",rownames(data.combined))
}
DefaultAssay(data.combined) <- "RNA"
clusters <- table(data.combined[['integrated_snn_res.0.3']])
clusters
for (i in 0:(length(clusters)-1)){
  print(i)
  custom.GO(markers03[[i+1]]$Gene,i,"res03",rownames(data.combined))
}

write.table(x = rownames(data.combined),file = "bg_genes.txt",quote = F,col.names = F,row.names = F)

custom.GO(markers03[[4+1]]$Gene,i,"res03",rownames(data.combined))
custom.GO(markers03[[5+1]]$Gene,i,"res03",rownames(data.combined))
markers03[[1]]
pdf(file = "./plots/umap_res03.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1 <- DimPlot(data.combined, reduction = "umap",
              group.by = "integrated_snn_res.0.3") + ggtitle("Louvain clusters")
p1
dev.off()
