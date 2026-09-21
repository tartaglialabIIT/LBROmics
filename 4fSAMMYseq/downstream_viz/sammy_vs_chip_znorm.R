# SAMMY-seq (S2S vs S3) vs ChIP-seq track figure + correlation heatmaps
# Outputs:
#  - Figure of example region on chr4:20,000,000-END with SAMMY mean +/- sd and ChIP tracks
#  - Heatmap of genome-wide Spearman correlation between SAMMY replicates and ChIP marks
#  - Heatmap of Spearman correlation among SAMMY replicates
#
# Requirements: Bioconductor packages Gviz, rtracklayer, GenomicRanges, IRanges, S4Vectors,
#               pheatmap, data.table, ggplot2
#
# Edit file paths below and run.

suppressPackageStartupMessages({
  library(Gviz)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(pheatmap)
  library(data.table)
  library(ggplot2)
})

# ---------------------------
# User inputs: change paths
# ---------------------------
wd <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/"
setwd(wd)

out_prefix <- "SAMMY_ChIP_figure"

cytodf <- read.table("./cytoBand.txt", header = FALSE, sep = "\t")

# ChIP bigwigs (from GSE96107) -- provide full paths
chip_files <- list(
  H3K9me3 = "./bonev_bigWig/GSE96107_NPC_H3K9me3.bw",
  H3K4me3 = "./bonev_bigWig/GSE96107_NPC_H3K4me3.bw",
  H3K27ac = "./bonev_bigWig/GSE96107_NPC_H3K27ac.bw"
)

# SAMMY-seq comparison bigwigs (three WT replicates)
sammy_files <- list(
  WT1 = "./comparisons/spp_mle/NPCwt_n1_S2SvsS3_mle.bigWig",
  WT2 = "./comparisons/spp_mle/NPCwt_n2_S2SvsS3_mle.bigWig",
  WT3 = "./comparisons/spp_mle/NPCwt_IRE_R3_S2SvsS3_mle.bigWig"
)

# Which chromosome/region to display for panel a
display_chr <- "chr6"
display_start <- 20000000
# We'll determine end from bigwig seqlengths; fallback to an arbitrary end if unavailable

# bin/window size for both display and genome-wide binning (match the methods: windowSize=1500)
bin_size <- 1500

# chromosomes to use for genome-wide correlation (autosomes + X if present)
# we'll fetch from one of the bigwigs

# ---------------------------
# Helper functions
# ---------------------------
# get seqinfo (seqlengths) from a bigWig file
get_seqinfo_from_bw <- function(bwfile) {
  bw <- BigWigFile(bwfile)
  si <- seqinfo(bw)
  return(si)
}

# create bins (GRanges) for a given region
make_bins_for_region <- function(chr, start, end, width) {
  starts <- seq(from = start, to = end, by = width)
  ends <- pmin(starts + width - 1, end)
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
}

# compute binned mean signal for a bigWig file over a set of bins (GRanges)
binned_signal_from_bw <- function(bwfile, bins) {
  # import intervals overlapping bins; import will bring 'score' values
  gr <- import(bwfile, which = bins)
  if (length(gr) == 0) {
    warning("No signal in ", bwfile, " for the provided bins")
    return(rep(0, length(bins)))
  }
  # compute coverage Rle weighted by score for each seqname present in bins
  cov <- coverage(gr, weight = "score")
  seqn <- as.character(seqnames(bins)[1])
  if (! (seqn %in% names(cov)) ) {
    warning("Chromosome ", seqn, " not present in coverage for file: ", bwfile)
    return(rep(0, length(bins)))
  }
  rle <- cov[[seqn]]
  # create IRanges of bin coordinates
  br <- ranges(bins)
  views <- Views(rle, start(br), end(br))
  means <- viewMeans(views)
  as.numeric(means)
}

# compute genome-wide tiled bins from seqinfo and selected chromosomes
make_genome_bins <- function(seqinfo, tilewidth, chromosomes = NULL) {
  if (is.null(chromosomes)) chromosomes <- names(seqlengths(seqinfo))
  # only keep chromosomes with known lengths
  seqlens <- seqlengths(seqinfo)[chromosomes]
  seqlens <- seqlens[!is.na(seqlens) & seqlens > 0]
  grlist <- GRangesList(lapply(names(seqlens), function(chr) {
    chrlen <- seqlens[[chr]]
    starts <- seq(1, chrlen, by = tilewidth)
    ends <- pmin(starts + tilewidth - 1, chrlen)
    GRanges(seqnames = chr, ranges = IRanges(starts, ends))
  }))
  unlist(grlist, use.names = FALSE)
}

# ---------------------------
# Determine seqinfo and display region end
# ---------------------------
# pick one bigWig (first SAMMY) to get seqinfo
first_bw <- sammy_files[[1]]
if (!file.exists(first_bw)) stop("First SAMMY bigwig not found: ", first_bw)
si <- get_seqinfo_from_bw(first_bw)

# determine display_end based on seqlengths if available
if (display_chr %in% names(seqlengths(si)) && !is.na(seqlengths(si)[display_chr])) {
  display_end <- seqlengths(si)[display_chr]
} else {
  # fallback: use start + 5 Mb or a large number
  display_end <- display_start + 5e6
  warning("Seqlength for ", display_chr, " not found. Using fallback end = ", display_end)
}

# make display bins
display_bins <- make_bins_for_region(display_chr, display_start, display_end, bin_size)

# Sample groups
sample_group1 <- "WT_NPC"
sample_group2 <- "MUT_NPC"

# Experimental IDs (for replicates)
exp_group1 <- c("NPCwt_n1", "NPCwt_n2", "NPCwt_IRE_R3")
exp_group2 <- c("NPCmutR1_Ire", "NPCmutR2_Ire", "NPCmut_IRE_R3")

pr <- matrix(list(exp_group1, exp_group2), nrow = 2)
rownames(pr) <- c(sample_group1, sample_group2)
colnames(cytodf) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

combination <- 1
x <- pr[1, combination][[1]]  

new_selection <- rdata$quantile_normalized_bins
mcols(new_selection) <- mcols(new_selection)[x]
sammy_dt <- DataTrack(
  range = new_selection,                             
  chromosome = "6",
  name = "SAMMY-seq S2SvsS3",
  genome = "mm10",
  type = c("a", "confint"),
  alpha.confint = 0.5,
  #alpha = 0.7,
  ucscChromosomeNames = FALSE,
  lwd.mountain = 2,
  cex.bands = 0.3,
  ylim = c(-0.8,0.8),
  window = 900,
  baseline = c(0),
  col.baseline = c("black"),
  col.axis = "black",
  col.yticks = "black",
  baseline.lty = "dotted"
)

make_binned_DataTrack <- function(file, name,col,  bins) {
  print(c(file,name,col))
  sig <- binned_signal_from_bw(file, bins)
  # Z-normalize ChIP tracks only
  if (name %in% names(chip_files)) {
    # robust center & scale
    med <- median(sig, na.rm = TRUE)
    madv <- mad(sig, constant = 1, na.rm = TRUE)  # raw MAD (no 1.4826)
    
    if (madv == 0 || is.na(madv)) {
      warning("Track ", chip, " has zero MAD; skipping robust z-normalization")
      next
    }
    
    # robust z-score
    sig <- (sig - med) / madv
    
    # clip negative values
    sig[sig < 0] <- 0
    
#    mu <- mean(sig, na.rm = TRUE)
#    sdv <- sd(sig, na.rm = TRUE)
#    sig <- (sig - mu) / sdv
  }
  gr  <- bins
  mcols(gr)$score <- sig
  DataTrack(
    range = gr,
    name = name,
    genome = "mm10",
    type = "histogram",
    col  = col,  
    lwd = 1.2,
    ylim = c(min(sig, na.rm = TRUE), max(sig, na.rm = TRUE)),
    add.yaxis = TRUE,
    ylab = name,
    col.axis = "black",
    col.yticks = "black"
  )
}

chip_bin_size <- 5000
chip_bins <- make_bins_for_region(display_chr, display_start, display_end, chip_bin_size)
chip_colors <- c(
  H3K27ac = "#C586E8",
  H3K4me3 = "#E85A2C",
  H3K9me3 = "#7C8C8E"
)

names(chip_files)

chip_tracks <- lapply(names(chip_files), function(nm) {
  make_binned_DataTrack(
    file = chip_files[[nm]],
    name = nm,
    col = chip_colors[[nm]],
    bins = chip_bins
  )
})
names(chip_tracks) <- names(chip_files)

# Chromosome ideogram and axis
ideo <- IdeogramTrack(chromosome = display_chr, genome = "mm10")
gaxis <- GenomeAxisTrack()

# Arrange the plot: ideogram, axis, SAMMY track and ChIP tracks stacked
chip_tracks
all_tracks <- c(list(ideo, gaxis,chip_tracks$H3K27ac,chip_tracks$H3K4me3, sammy_dt,chip_tracks$H3K9me3))#, chip_tracks)

#all_tracks <- c(list(ideo, gaxis, sammy_dt), chip_tracks)

pdf(paste0(out_prefix, "_panelA_chr6_regionZNORM.pdf"), width = 16, height = 6)
plotTracks(all_tracks, from = display_start, to = display_end, chromosome = display_chr,showId = TRUE,
           fontcolor.title = "black",
           col.sampleNames = "black",
           windowSize = bin_size, col.frame = "white", background.title = "transparent",add.yaxis = TRUE,
           fontsize = 12, legend = TRUE, ucscChromosomeNames = FALSE)
dev.off()
message("Saved panel A to ", paste0(out_prefix, "_panelA_chr4_region.pdf"))

# ---------------------------
# Genome-wide binning and correlation
# ---------------------------
# Determine genome-wide chromosomes to use (autosomes + X if present)
available_chrs <- names(seqlengths(si))
# prefer chr1-chr19, chrX for mouse (mm10) or human style chr1..chr22
preferred <- c(paste0("chr", 1:19), "chrX")
use_chrs <- intersect(preferred, available_chrs)
if (length(use_chrs) == 0) use_chrs <- available_chrs

message("Using chromosomes for genome-wide correlation: ", paste(use_chrs, collapse = ", "))
bin_size
genome_bins <- make_genome_bins(si, tilewidth = 50000, chromosomes = use_chrs)
message("Number of genome-wide bins: ", length(genome_bins))

# Function to compute binned signals for a named list of bigwig files
compute_binned_matrix <- function(bw_list, bins) {
  mat <- sapply(bw_list, function(f) {
    message("Binning ", f)
    binned_signal_from_bw(f, bins)
  })
  colnames(mat) <- names(bw_list)
  mat
}

# combine SAMMY (replicates) + ChIP (choose histone marks and lamin if available)
# user provided lamin? not in list — we'll proceed with H3 marks and H3K9me3 and Input
all_tracks_list <- c(sammy_files, chip_files[c("H3K27ac", "H3K4me3", "H3K9me3")])

# some of the files above might be duplicated; keep unique and named
all_tracks_list <- all_tracks_list[unique(names(all_tracks_list))]

# compute matrix (bins x tracks)
binned_mat <- compute_binned_matrix(all_tracks_list, genome_bins)

chip_names <- names(chip_files)


for (chip in chip_names) {
  if (chip %in% colnames(binned_mat)) {
    
    v <- binned_mat[, chip]
    
    # robust center & scale
    med <- median(v, na.rm = TRUE)
    madv <- mad(v, constant = 1, na.rm = TRUE)  # raw MAD (no 1.4826)
    
    if (madv == 0 || is.na(madv)) {
      warning("Track ", chip, " has zero MAD; skipping robust z-normalization")
      next
    }
    
    # robust z-score
    z <- (v - med) / madv
    
    # clip negative values
    z[z < 0] <- 0
    
    binned_mat[, chip] <- z
  }
}

# Remove bins that are all zeros or all NAs to avoid spurious correlations
keep <- apply(binned_mat, 1, function(r) { any(!is.na(r) & r != 0) })
message("Keeping ", sum(keep), " / ", nrow(binned_mat), " bins for correlation (non-zero in at least one track)")
binned_mat_filt <- binned_mat[keep, , drop = FALSE]

# compute Spearman correlations (pairwise complete observations)
corr_all <- cor(binned_mat_filt, method = "spearman", use = "pairwise.complete.obs")

# plot heatmap of correlations
pdf(paste0(out_prefix, "_panelB_correlation_alltracksZNORM.pdf"), width = 7, height = 6)
pheatmap(corr_all, main = "Genome-wide Spearman correlation (SAMMY + ChIP)",
         clustering_method = "complete", border_color = NA)
dev.off()
message("Saved correlation heatmap to ", paste0(out_prefix, "_panelB_correlation_alltracks.pdf"))

# SAMMY-only correlation
sammy_mat_genome <- binned_mat_filt[, intersect(names(sammy_files), colnames(binned_mat_filt)), drop = FALSE]
if (ncol(sammy_mat_genome) >= 2) {
  corr_sammy <- cor(sammy_mat_genome, method = "spearman", use = "pairwise.complete.obs")
  pdf(paste0(out_prefix, "_panelB_correlation_sammy_replicates.pdf"), width = 5, height = 4)
  pheatmap(corr_sammy, main = "SAMMY replicates Spearman correlation", border_color = NA)
  dev.off()
  message("Saved SAMMY replicates correlation heatmap to ", paste0(out_prefix, "_panelB_correlation_sammy_replicates.pdf"))
} else {
  warning("Not enough SAMMY replicates found in binned matrix to compute replicate correlation.")
}

# ---------------------------
# Done
# ---------------------------
message("All done. Output files with prefix: ", out_prefix)



library(pheatmap)
library(viridis)

library(RColorBrewer)
library(pheatmap)

# Define desired order
row_order <- c("H3K27ac", "H3K4me3", "H3K9me3", "WT1", "WT2", "WT3")
col_order <- c("WT1", "WT2", "WT3")

# Subset and reorder the correlation matrix
corr_sub <- corr_all[row_order, col_order]

# Define a Nature-style colormap from gold to forestgreen (-1 to 1)
my_colors <- colorRampPalette(c("#c7ae7e", "white", "#4c7882"))(100)

# Plot heatmap
pdf(paste0(out_prefix, "_panelB_correlation_heatmapZNORM.pdf"), width = 6, height = 5)
p <- pheatmap(corr_sub,
         color = my_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA,
         main = "Spearman correlation (SAMMY + ChIP)",
         fontsize_row = 10,
         fontsize_col = 10,
         display_numbers = TRUE,
         number_color = "black",
         breaks = seq(-1, 1, length.out = 101))
print(p)
dev.off()
message("Saved ordered correlation heatmap to ", paste0(out_prefix, "_panelB_correlation_heatmap.pdf"))
