#!/usr/bin/env Rscript

# Bulk RNA-seq DE analysis (DESeq2) + QC heatmaps/PCA + selected gene heatmaps
# Inputs (from GEO):
#   - raw_counts_NPC.txt   (NPC / differentiated)
#   - raw_counts_mESC.txt  (mESC / undifferentiated)
#
# Expected: count tables with first column = gene_id (Ensembl), remaining columns = samples.
# This script renames/reorders samples to the convention:
#   a_1  KO NPC  rep1
#   a_2  KO NPC  rep2
#   a_3  KO NPC  rep3
#   a_4  KO mESC rep1
#   a_5  KO mESC rep2
#   a_6  KO mESC rep3
#   a_7  WT NPC  rep1
#   a_8  WT NPC  rep2
#   a_9  WT mESC rep1
#   a_10 WT NPC  rep3
#   a_11 WT mESC rep2
#   a_12 WT mESC rep3

suppressPackageStartupMessages({
  library(DESeq2)
  library(rtracklayer)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
})

# -----------------------------
# User-configurable paths
# -----------------------------
counts_npc_path  <- "raw_counts_NPC.txt"
counts_mesc_path <- "raw_counts_mESC.txt"

# Put these in your repo (e.g., ./ref/) or point to where they are.
gtf_path      <- "ref/Mus_musculus.GRCm38.98.gtf"
gene_pos_path <- "ref/Mus_musculus.GRCm38.98_gen_pos.txt"  # tab-delimited, no header: gene_name chr start end

# Output folder (relative)
outdir <- "analysis_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# -----------------------------
# Helpers
# -----------------------------
read_counts <- function(path) {
  stopifnot(file.exists(path))
  x <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(ncol(x) >= 2)
  gene_id <- x[[1]]
  x <- x[, -1, drop = FALSE]
  rownames(x) <- gene_id
  # Ensure integer-ish matrix for DESeq2
  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  x
}

# If your count tables already have the final names, this will keep them.
# Otherwise, it will rename by ORDER to the convention above.
standardize_sample_names_by_order <- function(count_mat) {
  target <- c(
    "a_1 KO NPC rep1",
    "a_2 KO NPC rep2",
    "a_3 KO NPC rep3",
    "a_4 KO mESC rep1",
    "a_5 KO mESC rep2",
    "a_6 KO mESC rep3",
    "a_7 WT NPC rep1",
    "a_8 WT NPC rep2",
    "a_9 WT mESC rep1",
    "a_10 WT NPC rep3",
    "a_11 WT mESC rep2",
    "a_12 WT mESC rep3"
  )
  
  # If already standardized, do nothing
  if (all(colnames(count_mat) %in% target) && length(colnames(count_mat)) == 12) {
    return(count_mat[, target, drop = FALSE])
  }
  
  # Otherwise, if 12 samples, rename by order (most robust for public repo)
  if (ncol(count_mat) == 12) {
    message("Renaming samples by column order -> standardized a_* labels.")
    colnames(count_mat) <- target
    return(count_mat)
  }
  
  # If not 12 columns, fail loudly (better than silent wrong mapping)
  stop(
    "Unexpected number of sample columns: ", ncol(count_mat),
    ". Expected 12. If you have a different layout, add an explicit name mapping."
  )
}

build_metadata_from_colnames <- function(sample_names) {
  # sample_names like: "a_7 WT NPC rep1"
  tokens <- strsplit(sample_names, " ")
  df <- do.call(rbind, lapply(tokens, function(tk) {
    data.frame(
      sample   = paste(tk, collapse = " "),
      a_id     = tk[[1]],
      genotype = tk[[2]],        # WT / KO
      celltype = tk[[3]],        # NPC / mESC
      rep      = tk[[4]],        # rep1/2/3
      stringsAsFactors = FALSE
    )
  }))
  
  df$status <- ifelse(df$celltype == "NPC", "Diff", "Und")
  df$cell_line <- ifelse(df$genotype == "WT", "WT", "NT-KO")
  df$cell_line <- factor(df$cell_line, levels = c("WT", "NT-KO"))
  rownames(df) <- df$sample
  df
}

# Add gene symbols given an Ensembl->gene_name mapping
attach_symbols <- function(gtf_df, mat, id_colname = "gene_id") {
  df <- data.frame(ID = rownames(mat), mat, check.names = FALSE)
  out <- merge(unique(gtf_df[, c("gene_id", "gene_name")]),
               df,
               by.x = "gene_id",
               by.y = "ID",
               all = FALSE)
  out
}

# -----------------------------
# Load reference annotation
# -----------------------------
stopifnot(file.exists(gtf_path))
gtf <- rtracklayer::import(gtf_path)
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df[, c("gene_id", "gene_name")]
gtf_df <- gtf_df[!duplicated(gtf_df), ]
write.csv(gtf_df, "id_name_conversion.csv", row.names = FALSE)

# -----------------------------
# Load counts (NPC + mESC)
# -----------------------------
counts_npc  <- read_counts(file.path("..", counts_npc_path))
counts_mesc <- read_counts(file.path("..", counts_mesc_path))

# Merge to a unified matrix (union of genes; missing -> 0)
all_genes <- union(rownames(counts_npc), rownames(counts_mesc))
counts_all <- matrix(0L, nrow = length(all_genes), ncol = ncol(counts_npc) + ncol(counts_mesc),
                     dimnames = list(all_genes, c(colnames(counts_npc), colnames(counts_mesc))))

counts_all[rownames(counts_npc),  colnames(counts_npc)]  <- counts_npc
counts_all[rownames(counts_mesc), colnames(counts_mesc)] <- counts_mesc

# Standardize names + order
counts_all <- standardize_sample_names_by_order(counts_all)

# Metadata from standardized names
metadata <- build_metadata_from_colnames(colnames(counts_all))
write.csv(metadata, "metadata.csv", row.names = FALSE)

# Split
metadata.und  <- metadata[metadata$status == "Und", , drop = FALSE]
metadata.diff <- metadata[metadata$status == "Diff", , drop = FALSE]
counts.und    <- counts_all[, rownames(metadata.und),  drop = FALSE]
counts.diff   <- counts_all[, rownames(metadata.diff), drop = FALSE]

# -----------------------------
# DESeq2
# -----------------------------
dds.und <- DESeqDataSetFromMatrix(
  countData = counts.und,
  colData   = metadata.und,
  design    = ~ cell_line
)
dds.diff <- DESeqDataSetFromMatrix(
  countData = counts.diff,
  colData   = metadata.diff,
  design    = ~ cell_line
)

dds.und$cell_line  <- relevel(dds.und$cell_line,  "WT")
dds.diff$cell_line <- relevel(dds.diff$cell_line, "WT")

message("Genes before filtering: und=", nrow(dds.und), " diff=", nrow(dds.diff))
dds.und  <- dds.und[rowSums(counts(dds.und))  > 10, ]
dds.diff <- dds.diff[rowSums(counts(dds.diff)) > 10, ]
message("Genes after filtering:  und=", nrow(dds.und), " diff=", nrow(dds.diff))

dds.und2  <- DESeq(dds.und)
dds.diff2 <- DESeq(dds.diff)

# -----------------------------
# Normalized + raw counts export
# -----------------------------
dir.create("norm_counts", showWarnings = FALSE)
dir.create("raw_counts", showWarnings = FALSE)

norm_counts.und  <- log2(counts(dds.und2,  normalized = TRUE) + 1)
norm_counts.diff <- log2(counts(dds.diff2, normalized = TRUE) + 1)

write.table(
  attach_symbols(gtf_df, norm_counts.und),
  "norm_counts/normalized_counts_und.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(
  attach_symbols(gtf_df, norm_counts.diff),
  "norm_counts/normalized_counts_diff.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

raw_counts.und  <- counts(dds.und2)
raw_counts.diff <- counts(dds.diff2)

write.table(
  attach_symbols(gtf_df, raw_counts.und),
  "raw_counts/raw_counts_und.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(
  attach_symbols(gtf_df, raw_counts.diff),
  "raw_counts/raw_counts_diff.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

# -----------------------------
# QC: sample distance heatmaps + PCA
# -----------------------------
dir.create("plots", showWarnings = FALSE)

vsd.und  <- vst(dds.und2)
vsd.diff <- vst(dds.diff2)

sampleDistMatrix.und  <- as.matrix(dist(t(assay(vsd.und))))
sampleDistMatrix.diff <- as.matrix(dist(t(assay(vsd.diff))))

pheatmap(sampleDistMatrix.und,  fontsize = 12, filename = "plots/sample_distance_heatmap_und.pdf")
pheatmap(sampleDistMatrix.diff, fontsize = 12, filename = "plots/sample_distance_heatmap_diff.pdf")

pdf("plots/PCA_und.pdf", width = 6, height = 6)
print(plotPCA(vsd.und, intgroup = "cell_line"))
dev.off()

pdf("plots/PCA_diff.pdf", width = 6, height = 6)
print(plotPCA(vsd.diff, intgroup = "cell_line"))
dev.off()

# -----------------------------
# Differential expression (shrunken LFC) + exports
# -----------------------------
dir.create("DiffExp", showWarnings = FALSE)

# NOTE: lfcShrink(type="apeglm") requires apeglm installed.
de_shrink.und <- lfcShrink(dds = dds.und2,  coef = "cell_line_NT.KO_vs_WT", type = "apeglm")
de_shrink.diff <- lfcShrink(dds = dds.diff2, coef = "cell_line_NT.KO_vs_WT", type = "apeglm")

de_symbols.und <- merge(unique(gtf_df[, c("gene_id", "gene_name")]),
                        data.frame(ID = rownames(de_shrink.und), de_shrink.und, check.names = FALSE),
                        by.x = "gene_id", by.y = "ID", all = FALSE)

de_symbols.diff <- merge(unique(gtf_df[, c("gene_id", "gene_name")]),
                         data.frame(ID = rownames(de_shrink.diff), de_shrink.diff, check.names = FALSE),
                         by.x = "gene_id", by.y = "ID", all = FALSE)

write.table(de_symbols.und,  "DiffExp/deseq2_results_und.txt",  quote = FALSE, sep = "\t", row.names = FALSE)
write.table(de_symbols.diff, "DiffExp/deseq2_results_diff.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Significant sets (same thresholds as your original script)
de_symbols.und.sig <- na.omit(subset(de_symbols.und,  abs(log2FoldChange) > 1 & padj < 0.01))
de_symbols.diff.sig <- na.omit(subset(de_symbols.diff, abs(log2FoldChange) > 1 & padj < 0.01))

# -----------------------------
# Heatmaps for key gene sets (escapee, XCI silenced, reactivated, markers, top DE)
# -----------------------------
# Map Ensembl rownames -> gene_name for vst objects
rownames(vsd.diff) <- gtf_df$gene_name[match(rownames(vsd.diff), gtf_df$gene_id)]
rownames(vsd.und)  <- gtf_df$gene_name[match(rownames(vsd.und),  gtf_df$gene_id)]

# Column order for plots (same logic as original)
col.to.select.diff <- c("a_7 WT NPC rep1", "a_8 WT NPC rep2", "a_10 WT NPC rep3",
                        "a_1 KO NPC rep1", "a_2 KO NPC rep2", "a_3 KO NPC rep3")
col.to.select.und  <- c("a_9 WT mESC rep1", "a_11 WT mESC rep2", "a_12 WT mESC rep3",
                        "a_4 KO mESC rep1", "a_5 KO mESC rep2", "a_6 KO mESC rep3")

# Annotations
df.diff <- data.frame(Condition = metadata.diff$cell_line)
rownames(df.diff) <- rownames(metadata.diff)
df.und <- data.frame(Condition = metadata.und$cell_line)
rownames(df.und) <- rownames(metadata.und)

df.diff <- df.diff[col.to.select.diff, , drop = FALSE]
df.und  <- df.und[col.to.select.und,  , drop = FALSE]

df.diff$Condition <- factor(df.diff$Condition, levels = c("WT", "NT-KO"), labels = c("WT", "Lbr NT-KO"))
df.und$Condition  <- factor(df.und$Condition,  levels = c("WT", "NT-KO"), labels = c("WT", "Lbr NT-KO"))

mycolors <- list(Condition = setNames(c("#00BFC4", "#F8766D"), levels(df.diff$Condition)))

breaksList <- seq(-1.7, 1.7, by = 0.1)
hm_cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

# Gene sets
escapee.young <- c("Xist","Nkap","Car5b","Kdm5c","Mid1","Shroom4","Eif2s3x","Ddx3x")
escapee.mouse <- c("Sts","Pir","Ap1s2","Car5b","Txlng","Trappc2","Eif2s3x","Ddx3x","Kdm6a","Smc1a","Jpx")
escapee <- unique(c(escapee.young, escapee.mouse))

reactivated.mouse <- c("Ctps2","Nhs","Rab9","Asb9","Rbbp7","1810030O07Rik","Med14","Usp9x","Ndp","Nudt10","Srpx2","Tmem164","Pls3")

NPC_markers  <- c("Pax3","Pax6","Sox1","Sox2","Otx2","Ascl1","Smarca4","Msi1","Msi2","Nes")
mESC_markers <- c("Nanog","Pou5f1","Sox2")

# Escapee heatmaps (diff drives row order; und matches it)
escapee_in_diff <- escapee[escapee %in% rownames(vsd.diff)]
p.diff.escapee <- pheatmap(
  assay(vsd.diff)[escapee_in_diff, col.to.select.diff, drop = FALSE],
  cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.diff, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.diff.escapee, file = "DiffExp/heatmap_DIFF_escapee.pdf", device = "pdf")

ordered_escapee <- escapee_in_diff[p.diff.escapee$tree_row$order]
p.und.escapee <- pheatmap(
  assay(vsd.und)[ordered_escapee[ordered_escapee %in% rownames(vsd.und)], col.to.select.und, drop = FALSE],
  cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.und, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.und.escapee, file = "DiffExp/heatmap_UND_escapee.pdf", device = "pdf")

# NPC markers (diff)
npc_in_diff <- NPC_markers[NPC_markers %in% rownames(vsd.diff)]
p.diff.NPC <- pheatmap(
  assay(vsd.diff)[npc_in_diff, col.to.select.diff, drop = FALSE],
  cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.diff, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.diff.NPC, file = "DiffExp/heatmap_DIFF_NPCmarkers.pdf", device = "pdf")

# mESC markers (und + diff)
mesc_in_und <- mESC_markers[mESC_markers %in% rownames(vsd.und)]
p.und.mESC <- pheatmap(
  assay(vsd.und)[mesc_in_und, col.to.select.und, drop = FALSE],
  cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.und, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.und.mESC, file = "DiffExp/heatmap_UND_mESCmarkers.pdf", device = "pdf")

mesc_in_diff <- mESC_markers[mESC_markers %in% rownames(vsd.diff)]
p.diff.mESC <- pheatmap(
  assay(vsd.diff)[mesc_in_diff, col.to.select.diff, drop = FALSE],
  cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.diff, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.diff.mESC, file = "DiffExp/heatmap_DIFF_mESCmarkers.pdf", device = "pdf")

# Top 50 DE genes (by padj)
top_n <- 50
mESCs.SIG <- head(de_symbols.und.sig[order(de_symbols.und.sig$padj), ], top_n)
NPC.SIG   <- head(de_symbols.diff.sig[order(de_symbols.diff.sig$padj), ], top_n)

p.diff.top <- pheatmap(
  assay(vsd.diff)[NPC.SIG$gene_name[NPC.SIG$gene_name %in% rownames(vsd.diff)], col.to.select.diff, drop = FALSE],
  cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.diff, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.diff.top, file = "plots/heatmap_DIFF_DE_top_50_genes.pdf", device = "pdf", width = 6, height = 18)

p.und.top <- pheatmap(
  assay(vsd.und)[mESCs.SIG$gene_name[mESCs.SIG$gene_name %in% rownames(vsd.und)], col.to.select.und, drop = FALSE],
  cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
  annotation_col = df.und, annotation_colors = mycolors,
  show_rownames = TRUE, cellwidth = 10, cellheight = 10,
  color = hm_cols, breaks = breaksList
)
ggsave(p.und.top, file = "plots/heatmap_UND_DE_top_50_genes.pdf", device = "pdf", width = 6, height = 18)

# -----------------------------
# Supplementary tables: add chromosome to significant DE lists
# -----------------------------
dir.create("Supplementary_Tables", showWarnings = FALSE)

if (file.exists(gene_pos_path)) {
  gene_pos <- read.delim(gene_pos_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(gene_pos) <- c("gene_name", "chromosome", "start", "end")
  
  de_symbols.diff.sig$chromosome <- gene_pos$chromosome[match(de_symbols.diff.sig$gene_name, gene_pos$gene_name)]
  de_symbols.und.sig$chromosome  <- gene_pos$chromosome[match(de_symbols.und.sig$gene_name,  gene_pos$gene_name)]
  
  write.csv(de_symbols.diff.sig, "Supplementary_Tables/DE_genes_NPCs_day5.csv", row.names = FALSE)
  write.csv(de_symbols.und.sig,  "Supplementary_Tables/DE_genes_mESCs.csv",     row.names = FALSE)
} else {
  message("gene_pos_path not found (", gene_pos_path, "). Skipping chromosome-annotated supplementary tables.")
}

message("Done. Outputs written under: ", normalizePath(getwd()))