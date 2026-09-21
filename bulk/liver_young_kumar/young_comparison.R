wd <- "/Users/jonathan/Desktop/IIT/Cerase_single_cell/bulk_RNA_seq/analysis_results_MAY_2023/"
setwd(wd)

library(DESeq2)
library(rtracklayer)
library(readxl)

# Read the gtf file to map ensembl IDs to gene names
gtf <- rtracklayer::import('/Users/jonathan/Desktop/IIT/Cerase_single_cell/bulk_RNA_seq/Mus_musculus.GRCm38.98.gtf')
gtf_df=as.data.frame(gtf)
gtf_df <- gtf_df[,c("gene_id","gene_name")]
gtf_df <- gtf_df[!duplicated(gtf_df), ]

# Load the results of the DE in my NPC dataset
de_symbols.diff <- read.table("./DiffExp/deseq2_results_diff.txt",sep = "\t",header = T)

de_symbols.diff.sig <- de_symbols.diff[abs(de_symbols.diff$log2FoldChange) > 1,]
de_symbols.diff.sig <- de_symbols.diff.sig[de_symbols.diff.sig$padj < 0.01,]
de_symbols.diff.sig <- na.omit(de_symbols.diff.sig)

# Load the results of DE in female mice from Young et al 
my_data <- read_excel("../YoungData/Supplementary_Data/42003_2021_1944_MOESM4_ESM (1).xlsx", sheet = "fKO vs fWT")
my_data[my_data$FDR<0.05 & my_data$logFC>2 ,]

# Repeat DE on Young data with DESeq2

# NT-KO1 female : 5_17s001516
# NT-KO2 female: 6_17s001517
# WT1 female: 7_17s001518
# WT2 female: 8_17s001519

# NT-KO1 male: 2_17s001513
# NT-KO2 male: 3_17s001514
# WT1 male: 1_17s001512
# WT2 male: 4_17s001515

# Read the counts
counts.Young <- read.csv("../YoungData/GSE165447_Youngetalcountstable.csv",row.names = 1)
# Read the metadata
metadata.young <- read.csv("../YoungData/Young_Metadata.txt",sep = '\t',row.names = 1)
counts.Young <- counts.Young[,row.names(metadata.young)]


metadata.female <- metadata.young[metadata.young$sex == 'female',]
counts.female <- counts.Young[,row.names(metadata.female)]

dds <- DESeqDataSetFromMatrix(countData = counts.female,
                                   colData = metadata.female,
                                   design = ~ cell_line)

metadata.female$cell_line <- factor(metadata.female$cell_line, levels = c("NT-KO","WT"))

# Number of genes before filtering:
nrow(dds)

# Filter
dds <- dds[rowSums(counts(dds)) > 10, ]


# Number of genes left after low-count filtering:
nrow(dds)

dds2 <- DESeq(dds)
dds2$cell_line

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts.young <- log2(counts(dds2, normalized = TRUE)+1)

# add the gene symbols
norm_counts_young_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(norm_counts.young), norm_counts.young), by=1, all=F)

raw.counts.young <- counts(dds2)

raw_counts_young_symbols <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(raw.counts.young), raw.counts.young), by=1, all=F)


vsd <- vst(dds2)

vsd$cell_line <- factor(vsd$cell_line, levels = c("NT-KO","WT"))

resultsNames(dds2)

# Differential expression
de_shrink.young <- lfcShrink(dds = dds2,
                            coef="cell_line_WT_vs_NT.KO",
                            type="apeglm")

head(de_shrink.young)

# add the more comprehensive gene symbols to de_shrink
de_symbols.young <- merge(unique(gtf_df[,1:2]), data.frame(ID=rownames(de_shrink.young), de_shrink.young), by=1, all=F)
de_symbols.young.sig <- de_symbols.young[abs(de_symbols.young$log2FoldChange) > 1,]
de_symbols.young.sig <- de_symbols.young.sig[de_symbols.young.sig$padj < 0.01,]
de_symbols.young.sig <- na.omit(de_symbols.young.sig)
