library(DESeq2)
library(rtracklayer)

# Working directory: folder containing count tables / DiffExp outputs
# (set LBROMICS_B3_ROOT or run from that folder)
root <- Sys.getenv("LBROMICS_B3_ROOT", unset = ".")
setwd(root)

# GTF for gene names (GRCm38 / mm10)
gtf_file <- Sys.getenv("LBROMICS_GTF", unset = "mm10.gtf")
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)[, c("gene_id", "gene_name")]
gtf_df <- gtf_df[!duplicated(gtf_df), ]

# List count files
count_folder <- Sys.getenv("LBROMICS_B3_COUNTS", unset = "bulk_rnaseq_count_matrices")
ff <- list.files(path = count_folder, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(ff, read.table, skip = 4)

# Combine counts into a single matrix (use 4th column = reads for each gene)
counts <- as.data.frame(sapply(counts.files, function(x) x[, 4]))
rownames(counts) <- counts.files[[1]]$V1

# Extract sample names from filenames
sample_names <- gsub("_ReadsPerGene.out.tab$", "", basename(ff))
colnames(counts) <- sample_names

# Create metadata from filenames
metadata <- data.frame(
  sample = sample_names,
  cell_type = ifelse(grepl("ESC", sample_names), "ESC", "NPC"),
  genotype  = ifelse(grepl("WT", sample_names), "WT", "MUT")
)
rownames(metadata) <- metadata$sample

# Function to run DESeq2 for a subset
run_deseq <- function(count_matrix, meta_df, out_prefix) {
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = meta_df,
    design = ~ genotype
  )
  dds$genotype <- relevel(dds$genotype, "WT")
  
  # Filter low-count genes
  dds <- dds[rowSums(counts(dds)) > 10, ]
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Shrink log2 fold changes
  res <- lfcShrink(dds, coef="genotype_MUT_vs_WT", type="apeglm")
  
  # Add gene symbols
  res_symbols <- merge(gtf_df, data.frame(ID = rownames(res), res), by.x = "gene_id", by.y = "ID", all = FALSE)
  
  # Save results
  write.csv(res_symbols, paste0("./DiffExp/deseq2_results_", out_prefix, ".csv"), row.names = FALSE)
  
  dds
}

# Create output folder if it doesn't exist
dir.create("./DiffExp/", showWarnings = FALSE)

# Run DESeq2 separately for ESC and NPC
counts_ESC <- counts[, metadata$cell_type == "ESC"]
meta_ESC <- metadata[metadata$cell_type == "ESC", ]
dds_ESC <- run_deseq(counts_ESC, meta_ESC, "ESC")

counts_NPC <- counts[, metadata$cell_type == "NPC"]
meta_NPC <- metadata[metadata$cell_type == "NPC", ]
dds_NPC <- run_deseq(counts_NPC, meta_NPC, "NPC")

colnames(counts_NPC)


write.table(
  counts_ESC,
  file = "raw_counts_mESC.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


write.table(
  counts_NPC,
  file = "raw_counts_NPC.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(ggplotify)

# -----------------------------------------------------------
# Load results from CSV
# -----------------------------------------------------------
esc_res <- read.csv("./DiffExp/deseq2_results_ESC.csv")
npc_res <- read.csv("./DiffExp/deseq2_results_NPC.csv")

# Keep only unique gene names
esc_res <- esc_res[!duplicated(esc_res$gene_name), ]
npc_res <- npc_res[!duplicated(npc_res$gene_name), ]

# -----------------------------------------------------------
# Make PCA plots (ESC and NPC)
# -----------------------------------------------------------

make_pca_plot <- function(dds, prefix){
  
  vsd <- vst(dds, blind = TRUE)
  pca <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
  
  pca$genotype <- factor(pca$genotype,
                         levels = c("WT", "MUT"),
                         labels = c("WT", "Lbr NT-KO"))
  
  percentVar <- round(100 * attr(pca, "percentVar"))
  
  p <- ggplot(pca, aes(PC1, PC2, color = genotype)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("black", "red")) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste("PCA -", prefix)) +
    theme_bw(base_size = 16)
  
  ggsave(paste0("./DiffExp/PCA_", prefix, ".png"), p, width = 6, height = 5, dpi = 300)
  ggsave(paste0("./DiffExp/PCA_", prefix, ".pdf"), p, width = 6, height = 5)
  
  return(vsd)
}

vsd_ESC <- make_pca_plot(dds = dds_ESC, prefix = "ESC")
vsd_NPC <- make_pca_plot(dds = dds_NPC, prefix = "NPC")


# -----------------------------------------------------------
# Volcano plots
# -----------------------------------------------------------

make_volcano <- function(df, prefix){
  
  df$signif <- abs(df$log2FoldChange) > 1 & df$padj < 0.01
  
  p <- EnhancedVolcano(
    df,
    lab = df$gene_name,
    x = "log2FoldChange",
    y = "padj",
    title = paste("Volcano -", prefix),
    pCutoff = 0.01,
    FCcutoff = 1,
    labSize = 2.5,
    pointSize = 2,
    col = c("grey70", "grey70", "grey70", "red"),
    xlim = c(min(df$log2FoldChange, na.rm=TRUE), max(df$log2FoldChange, na.rm=TRUE))
  )
  
  ggsave(paste0("./DiffExp/Volcano_", prefix, ".png"), p, width = 7, height = 6, dpi = 300)
  ggsave(paste0("./DiffExp/Volcano_", prefix, ".pdf"), p, width = 7, height = 6)
  
}

make_volcano(esc_res, "ESC")
make_volcano(npc_res, "NPC")


# -----------------------------------------------------------
# Heatmaps: top 50 up + top 50 down
# -----------------------------------------------------------

make_heatmap <- function(dds, vsd, res_df, prefix){
  
  suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
  })
  
  # --------------------------------------------------
  # 1. Select significant DEGs
  # --------------------------------------------------
  sig <- res_df[res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1, ]
  
  # Top 50 up in MUT (log2FC > 0)
  top_mut_up <- head(sig[order(-sig$log2FoldChange), "gene_name"], 50)
  
  # Top 50 up in WT (log2FC < 0)
  top_wt_up  <- head(sig[order(sig$log2FoldChange), "gene_name"], 50)
  
  # Order: MUT-up first (top), WT-up next (bottom)
  genes_to_plot <- c(top_mut_up, top_wt_up)
  genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(vsd)]
  
  if (length(genes_to_plot) == 0) {
    message("No genes available for heatmap in ", prefix)
    return(NULL)
  }
  
  # --------------------------------------------------
  # 2. Force sample order: WT first, MUT second
  # --------------------------------------------------
  sample_info <- as.data.frame(colData(dds))
  sample_info$sample <- rownames(sample_info)
  
  sample_info <- sample_info[order(sample_info$genotype), ]  # WT first
  ordered_samples <- sample_info$sample
  
  # reorder VSD assay
  vsd_mat <- assay(vsd)[genes_to_plot, ordered_samples]
  
  # --------------------------------------------------
  # 3. Rename sample labels (ESC s1–s6, NPC s7–s12)
  # --------------------------------------------------
  if (prefix == "ESC") {
    new_labels <- paste0("s", 1:length(ordered_samples))
    mapping <- data.frame(original_sample = ordered_samples,
                          new_label = new_labels)
    write.csv(mapping, "./DiffExp/ESC_sample_mapping.csv", row.names = FALSE)
  } else if (prefix == "NPC") {
    new_labels <- paste0("s", (length(ordered_samples)+1):(2*length(ordered_samples)))
    mapping <- data.frame(original_sample = ordered_samples,
                          new_label = new_labels)
    write.csv(mapping, "./DiffExp/NPC_sample_mapping.csv", row.names = FALSE)
  } else {
    stop("prefix must be ESC or NPC")
  }
  
  colnames(vsd_mat) <- new_labels
  
  # --------------------------------------------------
  # 4. Build annotation (WT black, MUT red)
  # --------------------------------------------------
  annot <- data.frame(
    Condition = factor(
      sample_info$genotype,
      levels = c("WT", "MUT"),
      labels = c("WT", "Lbr NT-KO")
    )
  )
  rownames(annot) <- new_labels
  
  ann_colors <- list(
    Condition = c("WT" = "#000000", "Lbr NT-KO" = "#FF0000")
  )
  
  # --------------------------------------------------
  # 5. Heatmap colors
  # --------------------------------------------------
  breaksList <- seq(-2, 2, by = 0.1)
  hm_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(length(breaksList))
  
  # --------------------------------------------------
  # 6. Save PNG (no clipping)
  # --------------------------------------------------
  pheatmap(
    vsd_mat,
    cluster_rows = FALSE,      # prevent reordering of MUT-up / WT-up gene order
    cluster_cols = FALSE,      # keep WT → MUT order
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 10,
    color = hm_colors,
    breaks = breaksList,
    main = paste("Top DEGs -", prefix),
    filename = paste0("./DiffExp/Heatmap_", prefix, ".png"),
    width = 9,
    height = 12
  )
  
  # --------------------------------------------------
  # 7. Save PDF
  # --------------------------------------------------
  pheatmap(
    vsd_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 10,
    color = hm_colors,
    breaks = breaksList,
    main = paste("Top DEGs -", prefix),
    filename = paste0("./DiffExp/Heatmap_", prefix, ".pdf"),
    width = 9,
    height = 18
  )
  
  message("✔ Heatmap saved for ", prefix)
  message("✔ Mapping file saved for ", prefix)
}

make_heatmap(dds_ESC, vsd_ESC, esc_res, "ESC")
make_heatmap(dds_NPC, vsd_NPC, npc_res, "NPC")

ggsave(pesc, file="./DiffExp/PCA_ESC.pdf", device = "pdf")


# ============================================================
# HEATMAPS FOR NPCs - NEW MUTANT CLONE
#   1) Escapee genes
#   2) Neuronal progenitor markers
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(ggplotify)
})

dir.create("./DiffExp/", showWarnings = FALSE)

# ------------------------------------------------------------
# Make sure rownames of vsd_NPC are gene names
# ------------------------------------------------------------
gene_name_map <- gtf_df$gene_name[match(rownames(vsd_NPC), gtf_df$gene_id)]
rownames(vsd_NPC) <- gene_name_map

# remove rows without gene name
valid_rows <- !is.na(rownames(vsd_NPC)) & rownames(vsd_NPC) != ""
vsd_NPC <- vsd_NPC[valid_rows, ]

# if duplicated gene names exist, keep the first one
vsd_mat_full <- assay(vsd_NPC)
vsd_mat_full <- vsd_mat_full[!duplicated(rownames(vsd_mat_full)), ]

# ------------------------------------------------------------
# Order samples: WT first, then MUT
# ------------------------------------------------------------
sample_info <- as.data.frame(colData(dds_NPC))
sample_info$sample <- rownames(sample_info)

sample_info$genotype <- factor(sample_info$genotype, levels = c("WT", "MUT"))
sample_info <- sample_info[order(sample_info$genotype), ]

ordered_samples <- sample_info$sample
vsd_mat_full <- vsd_mat_full[, ordered_samples, drop = FALSE]

# ------------------------------------------------------------
# Rename columns if you want simple NPC labels
# ------------------------------------------------------------
new_labels <- paste0("s", 1:length(ordered_samples))
sample_mapping <- data.frame(
  original_sample = ordered_samples,
  new_label = new_labels
)
write.csv(sample_mapping, "./DiffExp/NPC_sample_mapping_for_special_heatmaps.csv", row.names = FALSE)

colnames(vsd_mat_full) <- new_labels

# ------------------------------------------------------------
# Column annotation
# ------------------------------------------------------------
annot <- data.frame(
  Condition = factor(
    sample_info$genotype,
    levels = c("WT", "MUT"),
    labels = c("WT", "Lbr NT-KO")
  )
)
rownames(annot) <- new_labels

ann_colors <- list(
  Condition = c("WT" = "#000000", "Lbr NT-KO" = "#FF0000")
)

# ------------------------------------------------------------
# Heatmap color scale
# ------------------------------------------------------------
breaksList <- seq(-2, 2, by = 0.1)
hm_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(length(breaksList))

# ============================================================
# 1) ESCAPEE GENES HEATMAP - NPC
# ============================================================

escapee.young <- c("Xist", "Nkap", "Car5b", "Kdm5c", "Mid1", "Shroom4", "Eif2s3x", "Ddx3x")
escapee.mouse <- c("Sts", "Pir", "Ap1s2", "Car5b", "Txlng", "Trappc2", "Eif2s3x", "Ddx3x", "Kdm6a", "Smc1a", "Jpx")
escapee <- unique(c(escapee.young, escapee.mouse))

escapee_present <- escapee[escapee %in% rownames(vsd_mat_full)]

cat("Escapee genes found in NPC matrix:\n")
print(escapee_present)

if (length(escapee_present) > 0) {
  
  # optionally restrict to genes present in DE results
  # escapee_present <- escapee_present[escapee_present %in% npc_res$gene_name]
  
  p.escapee <- pheatmap(
    vsd_mat_full[escapee_present, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 12,
    color = hm_colors,
    breaks = breaksList,
    main = "NPC - Escapee genes",
    filename = "./DiffExp/heatmap_NPC_escapee.pdf",
    width = 7,
    height = max(4, 0.35 * length(escapee_present) + 2)
  )
  
  pheatmap(
    vsd_mat_full[escapee_present, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 12,
    color = hm_colors,
    breaks = breaksList,
    main = "NPC - Escapee genes",
    filename = "./DiffExp/heatmap_NPC_escapee.png",
    width = 7,
    height = max(4, 0.35 * length(escapee_present) + 2)
  )
  
} else {
  message("No escapee genes found in vsd_NPC.")
}

# ============================================================
# 2) NEURONAL PROGENITOR MARKERS HEATMAP - NPC
# ============================================================

NPC_markers <- c("Pax3", "Pax6", "Sox1", "Sox2", "Otx2", "Ascl1", "Smarca4", "Msi1", "Msi2", "Nes")

npc_markers_present <- NPC_markers[NPC_markers %in% rownames(vsd_mat_full)]

cat("NPC marker genes found in NPC matrix:\n")
print(npc_markers_present)

if (length(npc_markers_present) > 0) {
  
  p.npc <- pheatmap(
    vsd_mat_full[npc_markers_present, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 12,
    color = hm_colors,
    breaks = breaksList,
    main = "NPC - Neuronal progenitor markers",
    filename = "./DiffExp/heatmap_NPC_neuronal_progenitor_markers.pdf",
    width = 7,
    height = max(4, 0.35 * length(npc_markers_present) + 2)
  )
  
  pheatmap(
    vsd_mat_full[npc_markers_present, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 12,
    color = hm_colors,
    breaks = breaksList,
    main = "NPC - Neuronal progenitor markers",
    filename = "./DiffExp/heatmap_NPC_neuronal_progenitor_markers.png",
    width = 7,
    height = max(4, 0.35 * length(npc_markers_present) + 2)
  )
  
} else {
  message("No NPC marker genes found in vsd_NPC.")
}

# ------------------------------------------------------------
# Optional: save which genes were actually plotted
# ------------------------------------------------------------
write.table(
  escapee_present,
  file = "./DiffExp/genes_plotted_escapee_NPC.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  npc_markers_present,
  file = "./DiffExp/genes_plotted_NPC_markers.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Done: NPC escapee and NPC marker heatmaps saved in ./DiffExp/")

npc.UP <- npc_res[(npc_res$log2FoldChange > 1 & npc_res$padj < 0.01), ]
npc.DOWN <- npc_res[(npc_res$log2FoldChange < -1 & npc_res$padj < 0.01), ]



# ============================================================
# 3) TOP 15 XCI-SILENCED DEGs HEATMAP - NPC
# ============================================================

# Read gene positions if not already loaded
gene_pos <- read.csv(
  Sys.getenv("LBROMICS_GENE_POS", unset = "../ref/Mus_musculus.GRCm38.98_gen_pos.txt"),
  sep = "\t",
  header = FALSE
)

colnames(gene_pos) <- c("gene_name", "chromosome", "start", "end")

# X chromosome genes
x_chr_genes <- gene_pos$gene_name[gene_pos$chromosome == "X"]

# Escapees already defined above:
# escapee.young
# escapee.mouse
# escapee

# XCI-silenced genes = X chromosome genes excluding escapees
xci_silenced_genes <- setdiff(x_chr_genes, escapee)

# Significant NPC DEGs
npc_sig <- npc_res[
  !is.na(npc_res$padj) &
    npc_res$padj < 0.01 &
    abs(npc_res$log2FoldChange) > 1,
]

# Keep only XCI-silenced genes among significant DEGs
npc_xcis_sig <- npc_sig[npc_sig$gene_name %in% xci_silenced_genes, ]

# Sort by adjusted p-value
npc_xcis_sig <- npc_xcis_sig[order(npc_xcis_sig$padj), ]

# Top 15
top15_xcis <- head(npc_xcis_sig$gene_name, 15)
top15_xcis <- top15_xcis[top15_xcis %in% rownames(vsd_mat_full)]

cat("Top 15 XCI-silenced DEGs found in NPC matrix:\n")
print(top15_xcis)

if (length(top15_xcis) > 0) {
  
  p.xcis <- pheatmap(
    vsd_mat_full[top15_xcis, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 12,
    color = hm_colors,
    breaks = breaksList,
    main = "NPC - Top 15 XCI-silenced DEGs",
    filename = "./DiffExp/heatmap_NPC_top15_XCI_silenced.pdf",
    width = 7,
    height = max(5, 0.35 * length(top15_xcis) + 2)
  )
  
  pheatmap(
    vsd_mat_full[top15_xcis, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    annotation_col = annot,
    annotation_colors = ann_colors,
    cellwidth = 14,
    cellheight = 12,
    color = hm_colors,
    breaks = breaksList,
    main = "NPC - Top 15 XCI-silenced DEGs",
    filename = "./DiffExp/heatmap_NPC_top15_XCI_silenced.png",
    width = 7,
    height = max(5, 0.35 * length(top15_xcis) + 2)
  )
  
} else {
  message("No XCI-silenced significant DEGs found in NPCs.")
}

# Optional: save the selected genes
write.table(
  top15_xcis,
  file = "./DiffExp/genes_plotted_top15_XCI_silenced_NPC.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



############### FIND DE genes on the X chromosome that should be highlighted in the karyoplot ##################
gene_pos <- read.csv(Sys.getenv("LBROMICS_GENE_POS", unset = "../ref/Mus_musculus.GRCm38.98_gen_pos.txt"),
                     sep= '\t',header=F)
gene.pos.X<- gene_pos[gene_pos$V2=="X",]

x.chr.genes <- gene.pos.X$V1

x.chr.genes.sig <- npc_sig[npc_sig$gene_name %in% x.chr.genes,]

gene.pos.X.sig <- gene.pos.X[gene.pos.X$V1 %in% x.chr.genes.sig$gene_name,]
colnames(gene.pos.X.sig) <- c("gene_name","chromosome","start","end")
gene.pos.X.sig$mid <- 0.5*(gene.pos.X.sig$end+gene.pos.X.sig$start)

gene.pos.X.sig$mid_Mb <- gene.pos.X.sig$mid/1000000
gene.pos.X.sig$bin <- cut(gene.pos.X.sig$mid_Mb,
                          breaks=hist(gene.pos.X.sig$mid_Mb,breaks = seq(0,175,5))$breaks)
gene.pos.X.sig$bin <- as.character(gene.pos.X.sig$bin)
library(dplyr)

bin.counts <- gene.pos.X.sig %>%
  group_by(bin) %>%
  summarise(n = n()) %>%
  as.data.frame()
bin.counts <- bin.counts[bin.counts$n >5,]
gene.pos.X.sig.karyo <- gene.pos.X.sig[gene.pos.X.sig$bin %in% bin.counts$bin,]

# Save the name of the genes to be highlighted
write.csv(x=gene.pos.X.sig.karyo,"genes_for_karyo_B3.csv")
