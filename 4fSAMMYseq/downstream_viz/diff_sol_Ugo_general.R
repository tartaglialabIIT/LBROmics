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

rdata$MUT_NPC_vs_WT_NPC_all_shifting_bins$
shifting.bins <- as.data.frame(rdata$MUT_NPC_vs_WT_NPC_all_shifting_bins)
shifting.bins$bintype

shifting.bins.X <- shifting.bins[shifting.bins$seqnames=="chrX",]
shifting.bins.X$bintype
write.csv(shifting.bins.X, "shiftingbinsX.csv", row.names=FALSE, quote=FALSE) 

# Load your quantile-normalized bins (GRanges)
bins <- rdata$quantile_normalized_bins



# Remove version numbers function
remove_version <- function(x) sub("\\..*$", "", x)

# Loop over each gene set and save to TXT
for (cat_name in names(rdata$genes)) {
  genes_cat <- rdata$genes[[cat_name]]
  genes_cat <- remove_version(genes_cat)  # remove version numbers
  file_name <- paste0(cat_name, ".txt")
  writeLines(genes_cat, file_name)
  message("Saved ", length(genes_cat), " genes to ", file_name)
}

# Cytoband file
cytodf <- read.table("../cytoBand.txt", header = FALSE, sep = "\t")

# Sample groups
sample_group1 <- "WT_NPC"
sample_group2 <- "MUT_NPC"

rdata$quantile_normalized_bins

# Experimental IDs (for replicates)
exp_group1 <- c("NPCwt_n1", "NPCwt_n2", "NPCwt_IRE_R3")
exp_group2 <- c("NPCmutR1_Ire", "NPCmutR2_Ire", "NPCmut_IRE_R3")

# Compared fractions (bin types)
bin_types_group1 <- c("S2S_up","S2S_down")
bin_types_group2 <- c("S3_up","S3_down")

# Combined bin types
all_bin_types <- c(bin_types_group1, bin_types_group2)

### =========================
### GENERAL DATA PROCESSING
### =========================

# Names in rdata
allgr_name <- paste0(sample_group2, "_allgr_", sample_group1)
genes_name  <- paste0(sample_group2, "_genes_", sample_group1)

# Load GRanges and genes list
all_gr <- rdata[[allgr_name]]
genes_list <- rdata[["genes"]]
counts <- lengths(genes_list)

# Create counts data frame
L2 <- sapply(names(counts), function(nm) {
  nm_upper <- toupper(nm)
  if (grepl("S2S", nm_upper) & grepl("UP", nm_upper)) "S2S_up" else
    if (grepl("S2S", nm_upper) & grepl("DOWN", nm_upper)) "S2S_down" else
      if (grepl("S3", nm_upper) & grepl("UP", nm_upper)) "S3_up" else "S3_down"
})

genes_count_df <- data.frame(
  L1 = factor(paste0(sample_group2, "_genes_", sample_group1), 
              levels = paste0(sample_group2, "_genes_", sample_group1)),
  L2 = factor(L2, levels = all_bin_types, ordered = TRUE),
  value = as.integer(counts)
)

# Ensure all categories are present
all_cats <- data.frame(
  L1 = paste0(sample_group2, "_genes_", sample_group1),
  L2 = all_bin_types
)
genes_count_df <- merge(all_cats, genes_count_df, all.x = TRUE)
genes_count_df$value[is.na(genes_count_df$value)] <- 0
genes_count_df$L2 <- factor(genes_count_df$L2, levels = all_bin_types, ordered = TRUE)

### =========================
### PLOT GENE COUNTS
### =========================

options(repr.plot.width =2, repr.plot.height = 6)
ggplot(genes_count_df, aes(x = L1, y = value, fill = L2)) +
  geom_col(width = 0.35, colour = "black") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 4, colour = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("S2S_up"="red","S2S_down"="orange","S3_up"="lightblue","S3_down"="darkblue"),
    breaks = all_bin_types
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())

### =========================
### PROCESS BIN DATA
### =========================
all_bin_types
paste0(sample_group2, "_", toupper("S2S_up"), "_", sample_group1)

rdata[["MUT_NPC_S2S_up_WT_NPC"]]

bin_dfs <- lapply(all_bin_types, function(bt) {
  df <- as.data.frame(rdata[[paste0(sample_group2, "_", bt, "_", sample_group1)]])
  df$bintype <- bt
  return(df)
})

allgr_df <- do.call(rbind, bin_dfs)
bins_df <- as.data.frame(rdata[["quantile_normalized_bins"]])

valid_btypes <- all_bin_types
long_all <- subset(allgr_df, bintype %in% valid_btypes, select = c("seqnames","start","end","width","bintype"))

chr_sizes <- aggregate(width ~ seqnames, bins_df, sum)
names(chr_sizes) <- c("seqnames","chr_tot_size")

df <- aggregate(width ~ seqnames + bintype, long_all, sum)
all_seq <- sort(unique(as.character(chr_sizes$seqnames)))
grid <- expand.grid(seqnames = all_seq, bintype = valid_btypes, stringsAsFactors = FALSE)
df <- merge(grid, df, by = c("seqnames","bintype"), all.x = TRUE)
df$width[is.na(df$width)] <- 0
df <- merge(df, chr_sizes, by = "seqnames", all.x = TRUE)
df$perc <- 100 * df$width / df$chr_tot_size
df$seqnames <- factor(df$seqnames, levels = c(paste0("chr",1:22),"chrX","chrY"), ordered = TRUE)
df$bintype  <- factor(df$bintype, levels = valid_btypes, ordered = TRUE)

### =========================
### PLOT SOLUBILITY COVERAGE
### =========================

options(repr.plot.width = 10, repr.plot.height = 6)
ggplot(df[df$seqnames != "chrY", ], aes(x = seqnames, y = perc, fill = bintype)) +
  geom_col() +
  coord_cartesian(ylim = c(0, 3)) +
  scale_y_continuous(breaks = seq(0, 3, 2)) +
  scale_fill_manual(values = c("red", "orange", "lightblue", "darkblue"), breaks = valid_btypes) +
  labs(
    #title = "Lbr NT-KO vs WT",
    x = NULL,
    y = "% of chromosome coverage",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 20),
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank()
  )

### =========================
### GVIZ TRACKS
### =========================

# Prepare comparison
pr <- matrix(list(exp_group1, exp_group2), nrow = 2)
rownames(pr) <- c(sample_group1, sample_group2)
colnames(cytodf) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

combination <- 1
x <- pr[1, combination][[1]]  
y <- pr[2, combination][[1]]  


xgroup <- sample_group1
ygroup <- sample_group2

contrast_col <- paste0(ygroup, "_vs_", xgroup)  
df1 <- as.data.frame(rdata[[combination]])  
stopifnot(all(c("seqnames","start","end", contrast_col) %in% colnames(df1)))
x1 <- makeGRangesFromDataFrame(df1[, c("seqnames","start","end", contrast_col)], keep.extra.columns = TRUE)

colnames(df1)
df1[, c("seqnames","start","end", contrast_col)]

df1X <- df1[df1$seqnames=="chrX",]

table(df1X$MUT_NPC_vs_WT_NPC)

### Function to create comparison track
func_plotting_comp_from_list <- function(z, xgroup, ygroup, xgrvect, ygrvect) {
  new_selection <- rdata$quantile_normalized_bins
  mcols(new_selection) <- mcols(new_selection)[z]
  dt_comparison <- DataTrack(
    range = new_selection,                             
    chromosome = "17",
    name = paste0(ygroup,"_vs_",xgroup),
    genome = "mm10",
    type = c("a", "confint"),
    alpha.confint = 0.5,
    alpha = 0.7,
    ucscChromosomeNames = FALSE,
    lwd.mountain = 2,
    cex.bands = 0.3,
    ylim = c(-0.8,0.8),
    window = 900,
    baseline = c(0),
    col.baseline = c("black"),
    baseline.lty = "dotted"
  )
  upcolors <- c(xgrvect, ygrvect)    
  displayPars(dt_comparison) <- list(
    groups = c(rep(paste0(1,ygroup), each = length(y)), rep(paste0(2,xgroup), each = length(x))),
    col = upcolors,
    legends = TRUE
  )
  return(dt_comparison)
}

# ORIGINAL FUNCTION WAS WITH GRADIENT
func_plottin_heatmap_from_gr <- function(gr_heat, ygroup, xgroup) {
  DataTrack(
    range = gr_heat,
    name  = paste0(xgroup, " vs ", ygroup),
    type  = "heatmap",
    # This is S3_down (-2), S3_up (-1), S2S_down (1), S2S_up (2)
    #col = mapped_colors,
    #gradient = c("darkblue","lightblue","orange","red"),
    gradient = c("lightblue","darkblue","red","orange"),
    showColorBar = TRUE,
    showAxis = FALSE,
    legend = TRUE
  )
  
}

table(x1$MUT_NPC_vs_WT_NPC)

comparis_plot1 <- func_plotting_comp_from_list(c(y, x), xgroup, ygroup, "red", "black")
comparis_heatmap1 <- func_plottin_heatmap_from_gr(x1, ygroup, xgroup)

comparis_heatmap1

options(repr.plot.width = 20, repr.plot.height = 4)
gtrack <- GenomeAxisTrack()

### =========================
### PLOT GVIZ TRACKS FOR ALL CHROMOSOMES
### =========================

for (i in c(seq(19),"X")) {
  ideoTrack <- IdeogramTrack(
    chromosome = paste0("chr", i),
    genome = "mm10",
    bands = cytodf,
    fontsize = 20
  )
  
  plotTracks(
    trackList = c(list(ideoTrack, gtrack), comparis_plot1, comparis_heatmap1),
    ucscChromosomeNames = FALSE,
    col.frame = "white",
    col.axis = "black",
    showId = TRUE,
    fontcolor.title = "black",
    col.sampleNames = "black",
    background.title = "transparent",
    windowSize = 1500,
    legend = TRUE
  )
}


