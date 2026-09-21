suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
  library(S4Vectors)
  library(ggplot2)
})

### =========================
### USER INPUTS (change for each comparison)
### =========================

# Path to working directory
wd <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/differential_solubility/"
setwd(wd)

# Input RDS file
rdata_file <- "rdata/S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds"
rdata <- readRDS(rdata_file)

my.genes <- rdata$genes

names(my.genes)


library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)

# ---- 1. Load DESeq2 results ----
res <- read_csv("/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/DESEQ2/DiffExp/deseq2_results_NPC.csv")

# Ensure gene names column is named "gene" (adjust if different)
colnames(res)[1] <- "gene"

# ---- 2. Define your named list of gene sets ----
gene_sets <- list(
  S2S_up   = my.genes$MUT_NPC_vs_WT_NPC_S2S_up,
  S2S_down = my.genes$MUT_NPC_vs_WT_NPC_S2S_down,
  S3_up    = my.genes$MUT_NPC_vs_WT_NPC_S3_up,
  S3_down  = my.genes$MUT_NPC_vs_WT_NPC_S3_down
)

# ---- 3. Extract log2FC values for each set ----
df_long <- imap_dfr(gene_sets, ~ {
  tibble(
    gene = .x,
    category = .y
  )
}) %>%
  left_join(res %>% select(gene, log2FoldChange), by = "gene") %>%
  drop_na(log2FoldChange)

# ---- 4. Boxplot with colors ----
p <- ggplot(df_long, aes(x = category, y = log2FoldChange, fill = category)) +
  geom_boxplot(outlier.size = 0.8) +
  theme_bw(base_size = 22) +
  scale_fill_manual(values = c(
    S2S_up   = "red",
    S2S_down = "orange",
    S3_up    = "lightblue",
    S3_down  = "darkblue"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log2 fold change") +
  xlab("") +
  ggtitle("Log2FC Lbr NT-KO vs WT from bulk RNA-seq of DSGs")

print(p)

# ---- 5. Statistical test: one-sample Wilcoxon test vs 0 ----
stats <- df_long %>%
  group_by(category) %>%
  summarize(
    p_value = wilcox.test(log2FoldChange, mu = 0)$p.value,
    median_log2FC = median(log2FoldChange)
  )

print(stats)
