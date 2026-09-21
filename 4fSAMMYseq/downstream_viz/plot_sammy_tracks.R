wd <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/differential_solubility/"
setwd(wd)

S2S_vsS3_WT_NPC_vs_MUT_NPC <- read_rds("./rdata/S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds")
names(S2S_vsS3_WT_NPC_vs_MUT_NPC)

all_shifting_bins <- S2S_vsS3_WT_NPC_vs_MUT_NPC$WT_NPC_vs_MUT_NPC_all_shifting_bins


S2S_vsS3_WT_NPC_vs_MUT_NPC$genes


library(Gviz)
library(rtracklayer)

bwfile <- "../NPCwt_n2_S2SvsS3_mle.bigWig"   # <- change this
chr <- "chr4"              # chromosome to plot
start <- 4000000           # start coordinate
end   <- 4080000           # end coordinate

# Load the bigWig signal as a track
bwTrack <- DataTrack(
  range = bwfile,
  genome = "mm10",       # or mm10, xenTro9, etc.
  chromosome = chr,
  type = "l",            # line plot
  name = "Signal"
)

# Genome axis track (x-axis)
axisTrack <- GenomeAxisTrack()

# Plot
plotTracks(list(axisTrack, bwTrack),
           from = start, to = end,
           col.line = "black")


library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

region <- GRanges("chr4", IRanges(30e6, 150e6))

bw <- import("../NPCwt_n2_S2SvsS3_mle.bigWig", which = region)

df <- as.data.frame(bw)

# bin into 1000 bins
nbins <- 1000
df$bin <- cut(df$start, breaks = nbins, labels = FALSE)
df_binned <- aggregate(score ~ bin, data = df, FUN = mean)

ggplot(df_binned, aes(x = bin, y = score)) +
  geom_line() +
  theme_bw() +
  labs(
    x = "chr4 (binned across 30–150 Mb)",
    y = "Avg Signal"
  )
