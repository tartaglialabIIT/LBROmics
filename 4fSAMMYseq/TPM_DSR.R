#!/usr/bin/env Rscript
# ============================================================
# TPM of genes overlapping DSRs (SAMMY-seq region-based method)
# ============================================================
#
# Usage:
#   Rscript TPM_DSR.R \
#     --rdata path/to/S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds \
#     --gtf path/to/mm10.gtf \
#     --tpm path/to/salmon.merged.gene_tpm.tsv \
#     --outdir ./FIGURES
#
# Or set LBROMICS_SAMMY_ROOT to a directory containing the relative paths
# used previously (differential_solubility/rdata/..., UCSC_mouse_genome/mm10.gtf,
# bulk_rnaseq_count_matrices/salmon_TPM/...).

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

root <- Sys.getenv("LBROMICS_SAMMY_ROOT", unset = getwd())
rdata_file <- get_arg("--rdata", file.path(root, "differential_solubility/rdata/S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds"))
gtf_file   <- get_arg("--gtf", file.path(root, "UCSC_mouse_genome/mm10.gtf"))
tpm_file   <- get_arg("--tpm", file.path(root, "bulk_rnaseq_count_matrices/salmon_TPM/salmon.merged.gene_tpm.tsv"))
outdir     <- get_arg("--outdir", file.path(root, "FIGURES"))

stopifnot(file.exists(rdata_file), file.exists(gtf_file), file.exists(tpm_file))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# Load SAMMY-seq DSRs
# -----------------------
rdata <- readRDS(rdata_file)

dsr_gr <- rdata$MUT_NPC_vs_WT_NPC_all_shifting_bins
table(dsr_gr$bintype)

length(rdata$genes$MUT_NPC_vs_WT_NPC_S2S_up)
length(rdata$genes$MUT_NPC_vs_WT_NPC_S2S_down)
length(rdata$genes$MUT_NPC_vs_WT_NPC_S3_up)
length(rdata$genes$MUT_NPC_vs_WT_NPC_S3_down)

# -----------------------
# Load gene annotation (UCSC GTF)
# -----------------------
gtf <- import(gtf_file)

# Keep only transcript features
tx_gr <- gtf[gtf$type == "transcript"]
tx_gr <- keepStandardChromosomes(tx_gr, pruning.mode = "coarse")

# Check required metadata
mcols(tx_gr)[, c("gene_id", "gene_name")] %>% head()

# -----------------------
# Build gene-level GRanges
# -----------------------
# -----------------------
# Build gene-level GRanges (UCSC GTF fix)
# -----------------------

# Split transcripts by gene
tx_by_gene <- split(tx_gr, tx_gr$gene_name)

# Reduce transcripts per gene
genes_grl <- reduce(tx_by_gene)

# Flatten to GRanges
genes_gr <- unlist(genes_grl, use.names = FALSE)

# Add gene_name from the GRangesList names
genes_gr$gene_name <- rep(
  names(genes_grl),
  elementNROWS(genes_grl)
)

# Optional: gene_id (first transcript per gene)
gene_id_map <- tapply(tx_gr$gene_id, tx_gr$gene_name, `[`, 1)
genes_gr$gene_id <- gene_id_map[genes_gr$gene_name]

# Sanity checks
class(genes_gr)
length(genes_gr)
head(genes_gr)

# -----------------------
# Overlap genes with DSRs
# -----------------------
hits <- findOverlaps(genes_gr, dsr_gr, ignore.strand = TRUE)
hits
gene_dsr_map <- data.frame(
  gene_name = genes_gr$gene_name[queryHits(hits)],
  gene_id   = genes_gr$gene_id[queryHits(hits)],
  bintype   = dsr_gr$bintype[subjectHits(hits)],
  stringsAsFactors = FALSE
)

# If a gene overlaps multiple DSRs, keep all (standard practice)
gene_dsr_map <- distinct(gene_dsr_map)

table(gene_dsr_map$bintype)

# -----------------------
# Load gene-level TPM matrix
# -----------------------
tpm_mat <- read.delim(tpm_file, stringsAsFactors = FALSE)

# -----------------------
# Subset TPMs to DSR-overlapping genes
# -----------------------
tpm_dsr <- tpm_mat %>%
  filter(gene_name %in% gene_dsr_map$gene_name) %>%
  left_join(gene_dsr_map, by = "gene_name")

# -----------------------
# Convert to long format
# -----------------------
tpm_long <- tpm_dsr %>%
  pivot_longer(
    cols = starts_with(c("ESC", "NPC")),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  mutate(
    CellType = ifelse(grepl("^ESC", Sample), "ESC", "NPC"),
    Condition = ifelse(grepl("[Mm][Uu][Tt]", Sample), "MUT", "WT"),
    log2TPM = log2(TPM + 1)
  )

# -----------------------
# Factor ordering
# -----------------------
tpm_long$bintype <- factor(
  tpm_long$bintype,
  levels = c("S2S_up", "S2S_down", "S3_down", "S3_up")
)

tpm_long$Group <- factor(
  paste0(tpm_long$bintype, "_", tpm_long$Condition),
  levels = c(
    "S2S_up_WT", "S2S_up_MUT",
    "S2S_down_WT", "S2S_down_MUT",
    "S3_down_WT", "S3_down_MUT",
    "S3_up_WT", "S3_up_MUT"
  )
)

# -----------------------
# Colors
# -----------------------
fill_colors <- c(
  "S2S_up_WT"    = "grey60",
  "S2S_up_MUT"  = "red",
  "S2S_down_WT" = "grey60",
  "S2S_down_MUT"= "orange",
  "S3_down_WT"  = "grey60",
  "S3_down_MUT" = "darkblue",
  "S3_up_WT"    = "grey60",
  "S3_up_MUT"   = "lightblue"
)

tpm_long <- tpm_long[tpm_long$CellType=="NPC",]

table(tpm_long$Condition)
table(tpm_long$Sample)
# -----------------------
# Boxplot
# -----------------------
out_pdf <- file.path(outdir, "new_plot_TPM_DSR.pdf")
pdf(out_pdf, height = 4, width = 7)
p_all <- ggplot(tpm_long, aes(x = Group, y = log2TPM, fill = Group)) +
  geom_boxplot(
    width = 0.65,
    color = "black",
    outlier.size = 0.6,
    outlier.alpha = 0.4
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(
    labels = c(
      "S2S\nUP", "", 
      "S2S\nDOWN", "", 
      "S3\nDOWN", "", 
      "S3\nUP", ""
    )
  )+
  scale_y_continuous(limits = c(-1, 15), expand = c(0, 0)) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  ylab("log2(TPM + 1)") +
  geom_signif(
    comparisons = list(
      c("S2S_up_WT", "S2S_up_MUT"),
      c("S2S_down_WT", "S2S_down_MUT"),
      c("S3_down_WT", "S3_down_MUT"),
      c("S3_up_WT", "S3_up_MUT")
    ),
    map_signif_level = TRUE,
    textsize = 5,
    tip_length = 0.01
  )

print(p_all)
dev.off()
message("Wrote ", normalizePath(out_pdf))
# -----------------------
# Wilcoxon tests
# -----------------------
wilcox_bycat <- tpm_long %>%
  group_by(CellType, bintype) %>%
  summarise(
    pval = wilcox.test(log2TPM ~ Condition)$p.value,
    n = n(),
    .groups = "drop"
  )

wilcox_bycat
