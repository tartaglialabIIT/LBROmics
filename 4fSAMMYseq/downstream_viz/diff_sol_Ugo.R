suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
  library(S4Vectors)
  library(ggplot2)
})

### REMEMBER TO CHANGE THE WD WITH YOUR PATH!

setwd("/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/differential_solubility/")
rdata <- readRDS("rdata/S2SvsS3_CM_CoVvsCM_NI_analysis.rds")
cytodf <- read.table("../cytoBand.txt", header = FALSE, sep = "\t")

### You can recover the hg38 cytodf in
### https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/

all_gr <- rdata[["CM_CoV_allgr_CM_NI"]]

serr_cols <- grep("_serrX2$", colnames(mcols(all_gr)), value = TRUE)

df <- as.data.frame(all_gr)
df$bin_id <- paste0(df$seqnames, "_", df$start, "_", df$end)
df$ratio  <- df[[serr_cols[1]]] / df[[serr_cols[2]]]

all_bins_df <- df[, c("seqnames", "start", "end", "bin_id", serr_cols, "ratio")]

all_gr <- rdata[[9]] #CM_CoV_allgr_CM_NI

CM_CoV_alldf_CM_NI <- as.data.frame(all_gr)

CM_CoV_alldf_CM_NI$bin_id <- paste0(
  CM_CoV_alldf_CM_NI$seqnames, "_",
  CM_CoV_alldf_CM_NI$start, "_",
  CM_CoV_alldf_CM_NI$end
)

genes_list <- rdata[["genes"]]
counts <- lengths(genes_list)

x1 <- setNames(data.frame(CM_CoV_genes_CM_NI = counts), "CM_CoV_genes_CM_NI")
x1






genes_list <- rdata[["genes"]]
counts <- lengths(genes_list)

L2 <- ifelse(grepl("S2S", names(counts)) & grepl("up", names(counts)), "S2S_up",
             ifelse(grepl("S2S", names(counts)) & grepl("down", names(counts)), "S2S_down",
                    ifelse(grepl("S3",  names(counts)) & grepl("up", names(counts)), "S3_up",
                           "S3_down")))

genes_count_df <- data.frame(
  L1    = factor("CM_CoV_genes_CM_NI", levels = "CM_CoV_genes_CM_NI"),
  L2    = factor(L2, levels = c("S2S_up","S2S_down","S3_up","S3_down"), ordered = TRUE),
  value = as.integer(counts),
  row.names = NULL
)

all_cats <- data.frame(L1 = "CM_CoV_genes_CM_NI",
                       L2 = c("S2S_up","S2S_down","S3_up","S3_down"))
genes_count_df <- merge(all_cats, genes_count_df, all.x = TRUE)
genes_count_df$value[is.na(genes_count_df$value)] <- 0

genes_count_df$L2 <- factor(genes_count_df$L2,
                            levels = c("S2S_up","S2S_down","S3_up","S3_down"),
                            ordered = TRUE)

options(repr.plot.width = 4, repr.plot.height = 6)

###Uncomment the line below if you want to generate the PDF
#pdf(file = "genes_counts.pdf", width = 4, height = 6)

ggplot(genes_count_df, aes(x = L1, y = value, fill = L2)) +
  geom_col(width = 0.35, colour = "black") +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5),
            size = 4, colour = "white", fontface = "bold") +
  scale_fill_manual(
    values = c(
      S2S_up   = "red",
      S2S_down = "orange",
      S3_up    = "lightblue",
      S3_down  = "darkblue"
    ),
    breaks = c("S2S_up","S2S_down","S3_up","S3_down")
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank()
  )




s2s_up   <- as.data.frame(rdata[["CM_CoV_S2S_up_CM_NI"]])
s2s_down <- as.data.frame(rdata[["CM_CoV_S2S_down_CM_NI"]])
s3_up    <- as.data.frame(rdata[["CM_CoV_S3_up_CM_NI"]])
s3_down  <- as.data.frame(rdata[["CM_CoV_S3_down_CM_NI"]])

s2s_up$bintype   <- "S2S_up"
s2s_down$bintype <- "S2S_down"
s3_up$bintype    <- "S3_up"
s3_down$bintype  <- "S3_down"

allgr_df <- rbind(s2s_up, s2s_down, s3_up, s3_down)

bins_df  <- as.data.frame(rdata[["quantile_normalized_bins"]])

valid_btypes <- c("S2S_up","S2S_down","S3_up","S3_down")
long_all <- subset(allgr_df, bintype %in% valid_btypes,
                   select = c("seqnames","start","end","width","bintype"))

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

options(repr.plot.width = 10, repr.plot.height = 6)

###Uncomment the line below if you want to generate the PDF
#pdf(file = "solubulty_coverage_per_chromosomes.pdf", width = 10, height = 6)

ggplot(df, aes(x = seqnames, y = perc, fill = bintype)) +
  geom_col() +
  coord_cartesian(ylim = c(0, 12)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_fill_manual(values = c("red","orange","lightblue","darkblue"),
                    breaks = valid_btypes) +
  labs(title = "CM_CoV_vs_CM_NI", x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank()
  )


pr <- matrix(list(
  c("CMR1","CMR2","CMR3"),
  c("CMR1CoV","CMR2CoV","CMR3CoV")
), nrow = 2)
rownames(pr) <- c("CM_NI", "CM_CoV")

colnames(cytodf) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

combination <- 1

x <- pr[1, combination][[1]]  
y <- pr[2, combination][[1]]  

xgroup <- "CM_NI"
ygroup <- "CM_CoV"

contrast_col <- paste0(ygroup, "_vs_", xgroup)  
df1 <- as.data.frame(rdata[[combination]])  

stopifnot(all(c("seqnames","start","end", contrast_col) %in% colnames(df1)))

x1 <- makeGRangesFromDataFrame(
  df1[, c("seqnames","start","end", contrast_col)],
  keep.extra.columns = TRUE
)

func_plotting_comp_from_list<- function(z,xgroup,ygroup,xgrvect,ygrvect){
  
  new_selection <-rdata$quantile_normalized_bins
  mcols(new_selection) <- mcols(new_selection)[z]
  dt_comparison<- DataTrack(range=new_selection,                             
                            chromosome="17",
                            name=paste0(ygroup,"_vs_",xgroup),
                            genome="hg38",
                            type=c("a", "confint"),
                            alpha.confint =0.5,
                            alpha=0.7,
                            ucscChromosomeNames=FALSE,
                            lwd.mountain = 2,
                            cex.bands=0.3,
                            ylim=c(-0.8,0.8),
                            window=900,baseline=c(0),col.baseline=c("black"),baseline.lty = "dotted")
  upcolors<-c(xgrvect,ygrvect)    
  displayPars(dt_comparison) <- list(
    groups =c(rep(paste0(1,ygroup),each = length(y) ),rep(paste0(2,xgroup),each = length(x) )),
    col= upcolors,
    legends=TRUE
  )
  return (dt_comparison)
}

func_plottin_heatmap<- function(combination,ygroup,xgroup){
  x<- makeGRangesFromDataFrame(as.data.frame(rdata[[combination]][[1]])[,c("seqnames","start","end",paste0(ygroup,"_vs_",xgroup))],keep.extra.columns = TRUE)
  
  prv_range_list <-  DataTrack(
    x,
    name = paste0(xgroup,"vs",ygroup),
    type = "heatmap",
    
    gradient = c( "darkblue","lightblue","orange","red" ),
    showColorBar = FALSE, showAxis = FALSE
  )
}

func_plottin_heatmap_from_gr <- function(gr_heat, ygroup, xgroup) {
  DataTrack(
    range = gr_heat,
    name  = paste0(xgroup, " vs ", ygroup),
    type  = "heatmap",
    gradient = c("darkblue","lightblue","orange","red"),
    showColorBar = TRUE, showAxis = FALSE, legend = TRUE
  )
}

comparis_plot1 <- func_plotting_comp_from_list(
  c(y, x),                 
  xgroup, ygroup,
  "#FDAC16", "#35AD78"
)
comparis_heatmap1 <- func_plottin_heatmap_from_gr(x1, ygroup, xgroup)

options(repr.plot.width=20, repr.plot.height=4)
gtrack <- GenomeAxisTrack()


### Uncomment the line below to save the plots to a PDF and keep dev.off() not commented
### Comment the PDF and comment ddev.off() to visualize the tracks on jupyter.
### Note: when pdf() is active, plots are written to the file and will not appear on screen until dev.off() is called.
#pdf(file = "S2SvsS3_consensus_gviz_chromosomes.pdf", width = 20, height = 4)

for (i in c(seq(22),"X","Y")) {          ### loop over all chromosomes
  # i <- 18                              ### uncomment to plot a specific chromosome
  
  ideoTrack <- IdeogramTrack(
    chromosome = paste0("chr", i),
    genome = "hg38",
    bands = cytodf,
    fontsize = 20
  )
  
  plotTracks(
    trackList = c(
      list(ideoTrack, gtrack),
      comparis_plot1,
      comparis_heatmap1
    ),
    ucscChromosomeNames = FALSE,
    # from = 35000000,                   ### uncomment to set start region
    # to = 75000000,                     ### uncomment to set end region
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

dev.off()