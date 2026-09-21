dir.create("/Users/jonathan/Desktop/IIT/Cerase_single_cell/bulk_RNA_seq/analysis_results_young_kumar/")
wd <- "/Users/jonathan/Desktop/IIT/Cerase_single_cell/bulk_RNA_seq/analysis_results_young_kumar/"
setwd(wd)

my.fisher <- function(DE.1,DE.2,total_genes){
  
  # Number of genes upregulated in both conditions
  overlap_genes <- length(intersect(DE.1, DE.2))
  
  # Number of genes upregulated only in list1
  only_list1 <- length(setdiff(DE.1, DE.2))
  
  # Number of genes upregulated only in list2
  only_list2 <- length(setdiff(DE.2, DE.1))
  
  # Number of genes not upregulated in either condition
  not_upregulated <- total_genes - (overlap_genes + only_list1 + only_list2)
  
  # Create a contingency table
  contingency_table <- matrix(c(overlap_genes, only_list1, only_list2, not_upregulated), nrow = 2)
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(contingency_table)
  
  fisher_result
}
  

plotPCA.mystyle <-  function(object, intgroup="condition", ntop=500, returnData=FALSE, pcs = c(1,2))
{
  
  stopifnot(length(pcs) == 2)    ### added this to check number of PCs ####
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  d <- data.frame(PC1=pca$x[,pcs[1]], PC2=pca$x[,pcs[2]], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcs]
    return(d)
  }
}

library(DESeq2)
library(rtracklayer)
library(readxl)

# Read the gtf file to map ensembl IDs to gene names
gtf <- rtracklayer::import('/Users/jonathan/Desktop/IIT/Cerase_single_cell/bulk_RNA_seq/Mus_musculus.GRCm38.98.gtf')
gtf_df=as.data.frame(gtf)
gtf_df <- gtf_df[,c("gene_id","gene_name")]
gtf_df <- gtf_df[!duplicated(gtf_df), ]
write.csv(gtf_df,"id_name_conversion.csv")

# Read the counts from Young et al (https://pubmed.ncbi.nlm.nih.gov/33846535/)
counts.Young <- read.csv("../YoungData/GSE165447_Youngetalcountstable.csv",row.names = 1)
# Read the metadata
metadata.young <- read.csv("../YoungData/Young_Metadata.txt",sep = '\t',row.names = 1)
counts.Young <- counts.Young[,row.names(metadata.young)]

# Create a column exp that we will use afterwards to control for batch effects in the 
# differential expression test
metadata.young$exp <- rep("young",length(row.names(metadata.young)))

# Read the counts from the kumar experiment
counts.kumar <- as.data.frame(read_excel("../Bulk_data_Kumar/countsNewTable.xlsx"))
rownames(counts.kumar) <- counts.kumar$Geneid
counts.kumar <- counts.kumar[ , !(colnames(counts.kumar) %in% c("Geneid"))]

# Read the metadata
metadata.kumar <- as.data.frame(read_excel("../Bulk_data_Kumar/RNA_analysis_LBR_KO_updated_version1.xlsx"))
rownames(metadata.kumar) <- metadata.kumar$FileNames
metadata.kumar <- metadata.kumar[ , !(colnames(metadata.kumar) %in% c("FileNames"))]

# Keep only the liver WT and Lbr NT-KO samples (4 samples)
metadata.kumar <- metadata.kumar[c("438","422","4612","4628"),]

counts.kumar <- counts.kumar[,row.names(metadata.kumar)]

# Create a column exp that we will use afterwards to control for batch effects in the 
# differential expression test
metadata.kumar$exp <- rep("kumar",length(row.names(metadata.kumar)))
metadata.kumar <- metadata.kumar[,c("Sex","Lbr_Condition","exp")]


metadata.kumar$Sex <- factor(metadata.kumar$Sex,levels = c("f","m"))
metadata.kumar$Lbr_Condition <- factor(metadata.kumar$Lbr_Condition,levels = c("WT","KO"))

library(plyr)

metadata.kumar$Sex <- revalue(metadata.kumar$Sex, c("f"="female", "m"="male"))
metadata.kumar$Lbr_Condition <- revalue(metadata.kumar$Lbr_Condition, c("WT"="WT", "KO"="NT-KO"))

colnames(metadata.kumar) <- colnames(metadata.young)
metadata.young$sex <- factor(metadata.young$sex) 
metadata.young$cell_line <- factor(metadata.young$cell_line) 

# Merge Young and Kumar data
full.counts <- merge(counts.Young, counts.kumar, by=0, all=TRUE) 
rownames(full.counts) <- full.counts$Row.names
full.counts <- full.counts[ , !(colnames(full.counts) %in% c("Row.names"))]

# Merge metadata
full.metadata <- rbind(metadata.young, metadata.kumar) 

all(colnames(full.counts)==rownames(full.metadata))

# Create a DESeq2 object
dds.full <- DESeqDataSetFromMatrix(countData = full.counts,
                                  colData = full.metadata,
                                  design = ~ exp + cell_line)
dds.full$cell_line <- relevel(dds.full$cell_line, "WT")

# Number of genes before filtering:
nrow(dds.full)

# Filter
dds.full <- dds.full[rowSums(counts(dds.full)) > 10, ]

# Number of genes left after low-count filtering:
nrow(dds.full)

dds.full2 <- DESeq(dds.full)

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts.full <- log2(counts(dds.full2, normalized = TRUE)+1)

# add the gene symbols
norm_counts_full_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(norm_counts.full), norm_counts.full), by=1, all=F)

# dir.create("./norm_counts/")

# write normalized counts to text file
# write.table(norm_counts_und_symbols, "./norm_counts/normalized_counts_und.txt", quote=F, col.names=T, row.names=F, sep="\t")

raw.counts.full <- counts(dds.full2)
raw_counts_full_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(raw.counts.full), raw.counts.full), by=1, all=F)

# dir.create("./raw_counts/")
# write.table(raw_counts_und_symbols, "./raw_counts/raw_counts_und.txt", quote=F, col.names=T, row.names=F, sep="\t")

vsd.full <- vst(dds.full2)

# load libraries pheatmap to create the heatmap plot
library(pheatmap)
library(ggplot2)
library(ggrepel)
# calculate between-sample distance matrix
sampleDistMatrix.full <- as.matrix(dist(t(assay(vsd.full))))

dir.create("./plots/")

pheatmap(sampleDistMatrix.full,cluster_cols = T,annotation = full.metadata,
         fontsize = 16, 
         filename = "./plots/sample_distance_heatmap_young_kumar.pdf",width = 8,height = 8)

# Make a PCA of the full dataset
pcaData <- plotPCA(vsd.full, intgroup=c("cell_line","sex","exp"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("./plots/PCA_young_kumar_all_12.pdf",width=6,height = 6)
p <- ggplot(pcaData, aes(x= PC1, y = PC2))+
  geom_point(size= 3, aes(shape=sex, fill=cell_line)) +
  scale_fill_manual(values = c('#00BFC4','#F8766D'),
                    guide = guide_legend(override.aes = list(shape = 22))) +
#  scale_fill_manual(values = c('#00BFC4','#F8766D'))+
  scale_shape_manual(values=c(21, 23, 25))+
  geom_text_repel(size= 3.5, aes(label=exp, colour=exp)) + 
  scale_color_manual(values = c("#E69F00", "#6756E9")) +
#  scale_color_discrete()+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Scaled data")
p
dev.off()

# Make a PCA of the full dataset
pcaData <- plotPCA.mystyle(vsd.full, intgroup=c("cell_line","sex","exp"), returnData=TRUE,pcs=c(1,4))

percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("./plots/PCA_young_kumar_all_14.pdf",width=6,height = 6)
p <- ggplot(pcaData, aes(x= PC1, y = PC2))+
  geom_point(size= 3, aes(shape=sex, fill=cell_line)) +
  scale_fill_manual(values = c('#00BFC4','#F8766D'),
                    guide = guide_legend(override.aes = list(shape = 22))) +
  #  scale_fill_manual(values = c('#00BFC4','#F8766D'))+
  scale_shape_manual(values=c(21, 23, 25))+
  geom_text_repel(size= 3.5, aes(label=exp, colour=exp)) + 
  scale_color_manual(values = c("#E69F00", "#6756E9")) +
  #  scale_color_discrete()+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC4: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Scaled data")
p
dev.off()


################ DIFFERENTIAL EXPRESSION ANALYSIS IN FEMALE SAMPLES ##############################

metadata.female <- full.metadata[full.metadata$sex == 'female',]
counts.female <- full.counts[,row.names(metadata.female)]

dds.female <- DESeqDataSetFromMatrix(countData = counts.female,
                              colData = metadata.female,
                              design = ~ exp + cell_line)

dds.female$cell_line <- relevel(dds.female$cell_line, "WT")

# Number of genes before filtering:
nrow(dds.female)

# Filter
dds.female <- dds.female[rowSums(counts(dds.female)) > 10, ]

# Number of genes left after low-count filtering:
nrow(dds.female)

dds2.female <- DESeq(dds.female)
dds2.female$cell_line

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts.female <- log2(counts(dds2.female, normalized = TRUE)+1)

# add the gene symbols
norm_counts_female_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(norm_counts.female), norm_counts.female), by=1, all=F)

raw.counts.female <- counts(dds2.female)

raw_counts_female_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(raw.counts.female), raw.counts.female), by=1, all=F)

vsd.female <- vst(dds2.female)

vsd.female$cell_line <- factor(vsd.female$cell_line, levels = c("NT-KO","WT"))

resultsNames(dds2.female)

# Differential expression
de_shrink.female <- lfcShrink(dds = dds2.female,
                           coef="cell_line_NT.KO_vs_WT",
                           type="apeglm")

head(de_shrink.female)

# add the more comprehensive gene symbols to de_shrink
de_symbols.female <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(de_shrink.female), de_shrink.female), by=1, all=F)

dir.create("./DiffExp/")

# write differential expression analysis result to a text file
write.table(de_symbols.female, "./DiffExp/deseq2_results_female_kumar_young.txt", quote=F, col.names=T, row.names=F, sep="\t")

de_symbols.female.sig <- de_symbols.female[abs(de_symbols.female$log2FoldChange) > 1,]
de_symbols.female.sig <- de_symbols.female.sig[de_symbols.female.sig$padj < 0.05,]
de_symbols.female.sig <- na.omit(de_symbols.female.sig)

length(rownames(de_symbols.female.sig))
length(rownames(de_symbols.female.sig[de_symbols.female.sig$log2FoldChange>0,]))
length(rownames(de_symbols.female.sig[de_symbols.female.sig$log2FoldChange<0,]))

# RESULTS DE in FEMALES from the full data (Young + Kumar):
# 4 genes up in Lbr NT-KO
# 4 genes down in Lbr NT-KO

# PCA of female samples
# Make a PCA of the full dataset
pcaData <- plotPCA(vsd.female, intgroup=c("cell_line","exp"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("./plots/PCA_young_kumar_female_12.pdf",width=6,height = 6)
p <- ggplot(pcaData, aes(x= PC1, y = PC2))+
  geom_point(size= 3, aes(shape=exp, fill=cell_line)) +
  scale_fill_manual(values = c('#00BFC4','#F8766D'),
                    guide = guide_legend(override.aes = list(shape = 22))) +
  #  scale_fill_manual(values = c('#00BFC4','#F8766D'))+
  scale_shape_manual(values=c(21, 23, 25))+
#  geom_text_repel(size= 3.5, aes(label=exp, colour=exp)) + 
#  scale_color_manual(values = c("#E69F00", "#6756E9")) +
  #  scale_color_discrete()+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Scaled data of female samples")
p
dev.off()

# Load the results of DE in female mice from Young et al 
young.female.DE <- as.data.frame(read_excel("../YoungData/Supplementary_Data/42003_2021_1944_MOESM4_ESM (1).xlsx",
                                            sheet = "fKO vs fWT"))
young.female.DE.sig <- young.female.DE[young.female.DE$FDR<0.05 & abs(young.female.DE$logFC)>1 ,]

length(rownames(young.female.DE.sig))
length(rownames(young.female.DE.sig[young.female.DE.sig$logFC>0,]))
length(rownames(young.female.DE.sig[young.female.DE.sig$logFC<0,]))

# RESULTS DE in FEMALES from Young data:
# 6 genes up in Lbr NT-KO
# 6 genes down in Lbr NT-KO

young.sig <- young.female.DE.sig$matchname
young.up <- young.female.DE.sig[young.female.DE.sig$logFC>0,]$matchname
young.down <- young.female.DE.sig[young.female.DE.sig$logFC<0,]$matchname

young.kumar.sig <- de_symbols.female.sig$gene_name
young.kumar.up <- de_symbols.female.sig[de_symbols.female.sig$log2FoldChange>0,]$gene_name
young.kumar.down <- de_symbols.female.sig[de_symbols.female.sig$log2FoldChange<0,]$gene_name

length(intersect(young.sig,young.kumar.sig))

library(VennDiagram)

g <- draw.pairwise.venn(area1=length(young.sig), area2=length(young.kumar.sig),
                   cross.area=length(intersect(young.sig,young.kumar.sig)), 
                   category=c("Young","Young&Kumar"))#,fill=c("Red","Yellow"))
pdf("./plots/venn_young_kumar_female.pdf",width = 6,height = 6)
plot.new()
grid.draw(g)
dev.off()

g.up <- draw.pairwise.venn(area1=length(young.up), area2=length(young.kumar.up),
                        cross.area=length(intersect(young.up,young.kumar.up)), 
                        category=c("Young","Young&Kumar"))#,fill=c("Red","Yellow"))
pdf("./plots/venn_young_kumar_female_up.pdf",width = 6,height = 6)
plot.new()
grid.draw(g.up)
dev.off()

g.down <- draw.pairwise.venn(area1=length(young.down), area2=length(young.kumar.down),
                        cross.area=length(intersect(young.down,young.kumar.down)), 
                        category=c("Young","Young&Kumar"))#,fill=c("Red","Yellow"))
pdf("./plots/venn_young_kumar_female_down.pdf",width = 6,height = 6)
plot.new()
grid.draw(g.down)
dev.off()

### HEATMAP FEMALE DATASET
df.female <- as.data.frame(colData(dds.female)[,c("cell_line")])
colnames(df.female) <- "Condition"
rownames(df.female) <- colnames(assay(vsd.female))
df.female$Condition <- factor(df.female$Condition,levels=c("WT","NT-KO"),labels=c("WT","Lbr NT-KO"))

mycolors <- c('#F8766D','#00BFC4')
names(mycolors) <- levels(df.female$Condition)
mycolors <- list(Condition = mycolors)

library(RColorBrewer)

pheatmap(assay(vsd.female),cluster_cols = T,
         scale = "column",
         show_rownames=F,
         annotation_col=df.female,
         # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         filename = "./plots/heatmap_young_kumar_female_scale_col.png",width = 6,height = 12)

### VOLCANO PLOT FEMALE SAMPLES
de_symbols.female$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de_symbols.female$diffexpressed[abs(de_symbols.female$log2FoldChange) > 1.0 & de_symbols.female$padj < 0.05] <- "DIFF"

max(de_symbols.female$log2FoldChange)
min(de_symbols.female$log2FoldChange)

de_symbols.female$logpadj <- -log10(de_symbols.female$padj)

# Re-plot but this time color the points with "diffexpressed"
pdf("./plots/volcano_plot_female_kumar_young2.pdf",width = 6,height = 6)

p <- ggplot(data=de_symbols.female, aes(x=log2FoldChange, y=logpadj, col=diffexpressed)) + geom_point() + theme_minimal()+
  theme(text = element_text(size = 18))

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c( "red", "black")
names(mycolors) <- c("DIFF", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3
dev.off()


################ DIFFERENTIAL EXPRESSION ANALYSIS IN MALE SAMPLES ##############################

metadata.male <- full.metadata[full.metadata$sex == 'male',]
counts.male <- full.counts[,row.names(metadata.male)]

dds.male <- DESeqDataSetFromMatrix(countData = counts.male,
                                   colData = metadata.male,
                                   design = ~ exp + cell_line)

dds.male$cell_line <- relevel(dds.male$cell_line, "WT")

# Number of genes before filtering:
nrow(dds.male)

# Filter
dds.male <- dds.male[rowSums(counts(dds.male)) > 10, ]

# Number of genes left after low-count filtering:
nrow(dds.male)

dds2.male <- DESeq(dds.male)
dds2.male$cell_line

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts.male <- log2(counts(dds2.male, normalized = TRUE)+1)

# add the gene symbols
norm_counts_male_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(norm_counts.male), norm_counts.male), by=1, all=F)

raw.counts.male <- counts(dds2.male)

raw_counts_male_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(raw.counts.male), raw.counts.male), by=1, all=F)

vsd.male <- vst(dds2.male)

vsd.male$cell_line <- factor(vsd.male$cell_line, levels = c("NT-KO","WT"))

resultsNames(dds2.male)

# Differential expression
de_shrink.male <- lfcShrink(dds = dds2.male,
                            coef="cell_line_NT.KO_vs_WT",
                            type="apeglm")

head(de_shrink.male)

# add the more comprehensive gene symbols to de_shrink
de_symbols.male <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(de_shrink.male), de_shrink.male), by=1, all=F)

dir.create("./DiffExp/")

# write differential expression analysis result to a text file
write.table(de_symbols.male, "./DiffExp/deseq2_results_male_kumar_young.txt", quote=F, col.names=T, row.names=F, sep="\t")

de_symbols.male.sig <- de_symbols.male[abs(de_symbols.male$log2FoldChange) > 1,]
de_symbols.male.sig <- de_symbols.male.sig[de_symbols.male.sig$padj < 0.05,]
de_symbols.male.sig <- na.omit(de_symbols.male.sig)

length(rownames(de_symbols.male.sig))
length(rownames(de_symbols.male.sig[de_symbols.male.sig$log2FoldChange>0,]))
length(rownames(de_symbols.male.sig[de_symbols.male.sig$log2FoldChange<0,]))

# RESULTS DE in MALES from the full data (Young + Kumar):
# 56 genes up in Lbr NT-KO
# 34 genes down in Lbr NT-KO

# PCA of male samples
# Make a PCA of the full dataset
pcaData <- plotPCA(vsd.male, intgroup=c("cell_line","exp"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("./plots/PCA_young_kumar_male_12.pdf",width=6,height = 6)
p <- ggplot(pcaData, aes(x= PC1, y = PC2))+
  geom_point(size= 3, aes(shape=exp, fill=cell_line)) +
  scale_fill_manual(values = c('#00BFC4','#F8766D'),
                    guide = guide_legend(override.aes = list(shape = 22))) +
  #  scale_fill_manual(values = c('#00BFC4','#F8766D'))+
  scale_shape_manual(values=c(21, 23, 25))+
  #  geom_text_repel(size= 3.5, aes(label=exp, colour=exp)) + 
  #  scale_color_manual(values = c("#E69F00", "#6756E9")) +
  #  scale_color_discrete()+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Scaled data of male samples")
p
dev.off()

# Load the results of DE in male mice from Young et al 
young.male.DE <- as.data.frame(read_excel("../YoungData/Supplementary_Data/42003_2021_1944_MOESM4_ESM (1).xlsx",
                                          sheet = "mKO vs mWT"))
young.male.DE.sig <- young.male.DE[young.male.DE$FDR<0.05 & abs(young.male.DE$logFC)>1 ,]

length(rownames(young.male.DE.sig))
length(rownames(young.male.DE.sig[young.male.DE.sig$logFC>0,]))
length(rownames(young.male.DE.sig[young.male.DE.sig$logFC<0,]))

# RESULTS DE in MALES from Young data:
# 111 genes up in Lbr NT-KO
# 98 genes down in Lbr NT-KO

young.sig <- young.male.DE.sig$matchname
young.up <- young.male.DE.sig[young.male.DE.sig$logFC>0,]$matchname
young.down <- young.male.DE.sig[young.male.DE.sig$logFC<0,]$matchname

young.kumar.sig <- de_symbols.male.sig$gene_name
young.kumar.up <- de_symbols.male.sig[de_symbols.male.sig$log2FoldChange>0,]$gene_name
young.kumar.down <- de_symbols.male.sig[de_symbols.male.sig$log2FoldChange<0,]$gene_name

length(intersect(young.sig,young.kumar.sig))

my.fisher(young.sig,young.kumar.sig,length(rownames(de_symbols.male)))
my.fisher(young.up,young.kumar.up,length(rownames(de_symbols.male[de_symbols.male$log2FoldChange>0,])))
my.fisher(young.down,young.kumar.down,length(rownames(de_symbols.male[de_symbols.male$log2FoldChange<0,])))


g <- draw.pairwise.venn(area1=length(young.sig), area2=length(young.kumar.sig),
                        cross.area=length(intersect(young.sig,young.kumar.sig)), 
                        category=c("Young","Young&Kumar"),fill = c("cornflowerblue", "darkorchid1"),
                        cex = 2,
                        cat.cex = 2,
                        cat.pos = c(275, 105),
                        cat.dist = 0.09,
                        cat.just = list(c(-1, -1), c(1, 1)))#,
#                        rotation.degree = 35)#,fill=c("Red","Yellow"))
pdf("./plots/venn_young_kumar_male.pdf",width = 6,height = 6)
plot.new()
title(main = "Male samples - All DEGs", cex.main = 2, line = -4, outer = TRUE)
grid.draw(g)
dev.off()

g.up <- draw.pairwise.venn(area1=length(young.up), area2=length(young.kumar.up),
                           cross.area=length(intersect(young.up,young.kumar.up)), 
                           category=c("Young","Young&Kumar"),
                           fill = c("cornflowerblue", "darkorchid1"),
                           cex = 2,
                           cat.cex = 2,
                           cat.pos = c(275, 105),
                           cat.dist = 0.09,
                           cat.just = list(c(-1, -1), c(1, 1)))
pdf("./plots/venn_young_kumar_male_up.pdf",width = 6,height = 6)
plot.new()
title(main = "Male samples - DEGs up in Lbr NT-KO", cex.main = 1.5, line = -4, outer = TRUE)
grid.draw(g.up)
dev.off()

g.down <- draw.pairwise.venn(area1=length(young.down), area2=length(young.kumar.down),
                             cross.area=length(intersect(young.down,young.kumar.down)), 
                             category=c("Young","Young&Kumar"),
                             fill = c("cornflowerblue", "darkorchid1"),
                             cex = 2,
                             cat.cex = 2,
                             cat.pos = c(275, 105),
                             cat.dist = 0.09,
                             cat.just = list(c(-1, -1), c(1, 1)))
pdf("./plots/venn_young_kumar_male_down.pdf",width = 6,height = 6)
plot.new()
title(main = "Male samples - DEGs down in Lbr NT-KO", cex.main = 1.5, line = -4, outer = TRUE)
grid.draw(g.down)
dev.off()

### HEATMAP MALE DATASET
df.male <- as.data.frame(colData(dds.male)[,c("cell_line")])
colnames(df.male) <- "Condition"
rownames(df.male) <- colnames(assay(vsd.male))
df.male$Condition <- factor(df.male$Condition,levels=c("WT","NT-KO"),labels=c("WT","Lbr NT-KO"))

mycolors <- c('#F8766D','#00BFC4')
names(mycolors) <- levels(df.male$Condition)
mycolors <- list(Condition = mycolors)

library(RColorBrewer)

length(rownames(assay(vsd.female)))
length(rownames(assay(vsd.male)))

pheatmap(assay(vsd.male),cluster_cols = T,
         scale = "column",
         show_rownames=F,
         annotation_col=df.male,
         # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         filename = "./plots/heatmap_young_kumar_male_scale_col.png",width = 6,height = 12)

### VOLCANO PLOT MALE SAMPLES
de_symbols.male$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de_symbols.male$diffexpressed[abs(de_symbols.male$log2FoldChange) > 1.0 & de_symbols.male$padj < 0.05] <- "DIFF"

# Clip values outside the ranges
de_symbols.male$log2FoldChange <- ramify::clip(de_symbols.male$log2FoldChange,-5,5)


de_symbols.male$logpadj <- -log10(de_symbols.male$padj)
de_symbols.male$logpadj <- ramify::clip(de_symbols.male$logpadj,0,10)
# Re-plot but this time color the points with "diffexpressed"
pdf("./plots/volcano_plot_male_kumar_young2.pdf",width = 6,height = 6)

p <- ggplot(data=de_symbols.male, aes(x=log2FoldChange, y=logpadj, col=diffexpressed)) + geom_point() + theme_minimal() +
  theme(text = element_text(size = 18))

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c( "red", "black")
names(mycolors) <- c("DIFF", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3
dev.off()

### SAVE DATA FOR GO MALES
write.table(x = de_symbols.male[de_symbols.male$diffexpressed=="DIFF" & de_symbols.male$log2FoldChange>0,]$gene_name,
            file = "./DiffExp/UP_male.txt",sep = "\t",quote = F,row.names = F)
write.table(x = de_symbols.male[de_symbols.male$diffexpressed=="DIFF" & de_symbols.male$log2FoldChange<0,]$gene_name,
            file = "./DiffExp/DOWN_male.txt",sep = "\t",quote = F,row.names = F)
write.table(x = de_symbols.male$gene_name,
            file = "./DiffExp/bg_male.txt",sep = "\t",quote = F,row.names = F)
