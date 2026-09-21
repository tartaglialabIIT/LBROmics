# --- Libraries ---
library(GenomicRanges)
library(rtracklayer)
library(preprocessCore)
library(Gviz)
library(dplyr)

# Set LBROMICS_SAMMY_ROOT or run from the nf-core comparisons output directory
sammy_root <- Sys.getenv("LBROMICS_SAMMY_ROOT", unset = ".")
setwd(file.path(sammy_root, "comparisons"))

# --- 1. Import bigWig files for NPC S2SvsS3 ---
npc_wt_files   <- c("./spp_mle/NPCwt_n1_S2SvsS3_mle.bigWig", "./spp_mle/NPCwt_n2_S2SvsS3_mle.bigWig","./spp_mle/NPCwt_IRE_R3_S2SvsS3_mle.bigWig")
npc_mut_files  <- c("./spp_mle/NPCmutR2_Ire_S2SvsS3_mle.bigWig", "./spp_mle/NPCmut_IRE_R3_S2SvsS3_mle.bigWig","./spp_mle/NPCmutR1_Ire_S2SvsS3_mle.bigWig")

# Function to import and resize as GRanges with 50kb bins
import_bigwig <- function(files, bin_size = 50000) {
  gr_list <- lapply(files, function(f) {
    gr <- import(f, format = "bigWig")
    # Resize bins to 50kb
    gr <- resize(gr, width = bin_size, fix = "start")
    gr
  })
  return(gr_list)
}

npc_wt_gr  <- import_bigwig(npc_wt_files)
npc_mut_gr <- import_bigwig(npc_mut_files)

# --- 2. Restrict to X chromosome ---
subset_chrX <- function(gr_list) lapply(gr_list, function(gr) gr[seqnames(gr) == "chrX"])
npc_wt_gr  <- subset_chrX(npc_wt_gr)
npc_mut_gr <- subset_chrX(npc_mut_gr)

# Combine all GRanges lists
all_gr_list <- c(npc_wt_gr, npc_mut_gr)

# Find the minimum length across ALL samples
min_len <- min(sapply(all_gr_list, length))

# Function to extract only min_len scores
gr_to_matrix <- function(gr_list, min_len) {
  mat <- sapply(gr_list, function(gr) gr$score[1:min_len])
  return(as.matrix(mat))
}

# Extract matrices with the same number of bins
wt_mat  <- gr_to_matrix(npc_wt_gr, min_len)
mut_mat <- gr_to_matrix(npc_mut_gr, min_len)

# --- 4. Quantile normalization ---
wt_mat_norm  <- normalize.quantiles(wt_mat)
mut_mat_norm <- normalize.quantiles(mut_mat)

# --- 5. Consensus track: mean + SE ---
wt_mean <- rowMeans(wt_mat_norm)
wt_se   <- apply(wt_mat_norm, 1, sd) / sqrt(ncol(wt_mat_norm))

mut_mean <- rowMeans(mut_mat_norm)
mut_se   <- apply(mut_mat_norm, 1, sd) / sqrt(ncol(mut_mat_norm))


# --- 6. GRanges as before ---
n_bins <- length(mut_mean)
chrX_bins <- npc_mut_gr[[1]][1:n_bins]

wt_consensus <- chrX_bins
wt_consensus$score <- wt_mean
wt_consensus$upper <- wt_mean + wt_se
wt_consensus$lower <- wt_mean - wt_se

mut_consensus <- chrX_bins
mut_consensus$score <- mut_mean
mut_consensus$upper <- mut_mean + mut_se
mut_consensus$lower <- mut_mean - mut_se


# --- 7. TRACKS ---

library(Gviz)

# --- 1. Basic tracks ---
ideo_track <- IdeogramTrack(genome = "mm10", chromosome = "chrX")
axis_track <- GenomeAxisTrack()

# --- 2. WT & MUT consensus tracks with CI ---
wt_track <- DataTrack(
  range = wt_consensus,
  data = "score",
  type = "a",                 # area plot
  name = "WT NPC",
  col = "black",
  fill = "lightgrey",
  window = "auto",
  confint = cbind(wt_consensus$lower, wt_consensus$upper),  # shaded CI,
  legend=T
)

mut_track <- DataTrack(
  range = mut_consensus,
  data = "score",
  type = "a",
  name = "Lbr NT-KO",
  col = "red",
  fill = "pink",
  window = "auto",
  confint = cbind(mut_consensus$lower, mut_consensus$upper),
  legend = T
)

# --- 3. Overlay tracks ---
overlay_track <- OverlayTrack(trackList = list(wt_track, mut_track))

# Plot tracks
ideo_track <- IdeogramTrack(genome = "mm10", chromosome = "chrX")
axis <- GenomeAxisTrack()

# WT CI polygon
wt_ci <- DataTrack(
  range = wt_consensus,
  data = wt_consensus$upper,
  baseline = wt_consensus$lower,
  type = "a",
  col = NA,
  fill = "lightblue",
  name = "WT NPC"
)
wt_line <- DataTrack(
  range = wt_consensus,
  data = "score",
  type = "l",
  col = "blue",
  lwd = 2,
  name = "WT NPC"
)

# MUT CI polygon
mut_ci <- DataTrack(
  range = mut_consensus,
  data = mut_consensus$upper,
  baseline = mut_consensus$lower,
  type = "a",
  col = NA,
  fill = "pink",
  name = "Lbr NT-KO"
)
mut_line <- DataTrack(
  range = mut_consensus,
  data = "score",
  type = "l",
  col = "red",
  lwd = 2,
  name = "Lbr NT-KO"
)

# Overlay WT + MUT
overlay <- OverlayTrack(trackList = list(wt_ci, mut_ci, wt_line, mut_line))

plotTracks(
  list(ideo_track, axis, overlay),
  chromosome = "chrX",
  from = start(chrX_bins)[1],
  to   = end(chrX_bins)[length(chrX_bins)],
  main = "SAMMY-seq consensus tracks - NPC S2SvsS3",
  legend = TRUE,
  legendPos = "topright"
)
