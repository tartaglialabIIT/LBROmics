library(tidyverse)
library(data.table)
library(ggplot2)
library(readr)

wd <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/"
setwd(wd)

sammy <- readRDS("./differential_solubility/rdata/S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds")
names(sammy)

sammy$MUT_NPC_vs_WT_NPC_all_shifting_bins

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(stringr)

### -------------------------
### 1. Load SAMMY GRanges
### -------------------------
gr <- sammy$MUT_NPC_vs_WT_NPC_all_shifting_bins

### -------------------------
### 2. Load GTF to map ranges → gene symbol
### -------------------------
gtf <- rtracklayer::import('/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/DESEQ2/mm10.gtf')

# keep only "gene" features with gene_id + gene_name
gtf_df <- gtf[gtf$type == "gene"]

gtf_annot <- data.frame(
  seqnames = as.character(seqnames(gtf_df)),
  start = start(gtf_df),
  end = end(gtf_df),
  gene_id = gtf_df$gene_id,
  gene_name = gtf_df$gene_name
)

### Overlap SAMMY genomic regions with GTF genes
hits <- findOverlaps(gr, gtf)
hits
sammy_df <- as.data.frame(gr)
class(sammy_df)
sammy_df$gene_name <- NA
sammy_df$gene_id <- NA

colnames(sammy_df)

sammy_df$gene_name[queryHits(hits)] <- gtf$gene_name[subjectHits(hits)]
sammy_df$gene_name

#sammy_df$gene_id[queryHits(hits)]   <- gtf_df$gene_id[subjectHits(hits)]

# keep unique by region (some bins may overlap >1 gene)
sammy_df <- sammy_df %>% distinct()

### Keep only bins that map to a gene
sammy_df <- sammy_df %>% filter(!is.na(gene_name))

### -------------------------
### 3. Load bulk DEGs
### -------------------------
bulk <- fread("./DESEQ2/DiffExp/deseq2_results_NPC.csv")

# Expecting columns: gene, log2FoldChange, padj
# Rename if needed
colnames(bulk) <- str_replace(colnames(bulk), "^Gene$", "gene_name")

### -------------------------
### 4. Merge SAMMY + bulk
### -------------------------
merged <- sammy_df %>%
  select(gene_name, cohen.estimate) %>%
  left_join(
    bulk %>% select(gene_name, log2FoldChange, padj),
    by = "gene_name"
  ) %>%
  filter(!is.na(log2FoldChange))

### -------------------------
### 5. Scatterplot
### -------------------------

ggplot(merged, aes(x = log2FoldChange, y = cohen.estimate)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 14) +
  labs(
    x = "Bulk log2 Fold Change",
    y = "SAMMY Cohen's d",
    title = "NPC - SAMMY Δsolubility (Cohen's d) vs Bulk log2FC (S2SvsS3)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black")
  )
