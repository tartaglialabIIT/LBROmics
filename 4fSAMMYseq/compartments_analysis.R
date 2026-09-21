### importing libraries
suppressPackageStartupMessages({
  require( data.table )  
  require( GenomicRanges )
  library(ggplot2) 
  require( factoextra )
})


## CONFIGURAZIONE
#bl_gr<- rtracklayer::import('/storage-daredevil/sammyseq_nfcore/Data/MEF/asset_mm10/ENCFF547MET.bed')
bed_dir <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/compartments/"
bedgraph_eigen_dir <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/compartments/"

## Importing bed compartments and bedgraphs:
  ### Defining sample groups
  cellular_types_list <- list(
#    MEFwt= c('NPCwt_n1','NPCwt_n2','NPCwt_IRE_R3'),
#    MEFACKO= c('NPCmutR1_Ire','NPCmutR2_Ire','NPCmut_IRE_R3')
    NPCwt = c('NPCwt_n1', 'NPCwt_n2','NPCwt_IRE_R3'),
    NPCmut = c('NPCmutR1_Ire', 'NPCmutR2_Ire', 'NPCmut_IRE_R3')
      )
### Importing bed compartments
list_bed<- list()
for (name in unlist(cellular_types_list)){
  print(name)
  # df_list<- lapply( 1:19, function( chr_n ){
  
  #     names<- paste0(name,"_","chr",chr_n,"_compartments_merged.tsv")
  #     df_chr<-read.table(paste0("",name,"_","chr",chr_n,"_compartments_merged.tsv"),header = T,sep = "\t")
  #     #
  #     #chr_gr<- makeGRangesFromDataFrame( df_chr,keep.extra.columns=T)
  #     #granges_list<- list(granges_list ,chr_gr)
  
  #     #return(chr_gr)
  #     return(df_chr)
  bed<- rtracklayer::import.bed(paste0(bed_dir,name,"_combined_compartments.bed"))
  if (min(start(bed)) == 3) {
    start(bed) <- start(bed) - 2
    end(bed) <- end(bed) - 1
  }
  #list_bed<- c(list_bed,bed)
  list_bed[[name]] <- bed  
}
names(list_bed)<- unlist(cellular_types_list)

df_bins<- as.data.frame(list_bed[['NPCwt_n1']])[c("seqnames","start","end")]   ### CAMBIARE NOME CON QUELLO DEL PRIMO INDICATO in cellular_types_list
#df_bins$start<-df_bins$start-2
for (name in names(list_bed)){
  df_new_sample<- as.data.frame(list_bed[[name]])[c("seqnames","start","end","name")]
  colnames(df_new_sample)<- c("seqnames","start","end",name)
  #df_new_sample$start<-df_new_sample$start-2
  df_bins<- merge(df_bins,df_new_sample,all = TRUE,by = c("seqnames","start","end"))
  
}

chrom_order <- paste0("chr", c(1:19,"X"))   ### CAMBIARE IN BASE AL GENOMA E AI CROMOSOMI CHE SI DEVE ANALIZZARE
df_bins$seqnames <- factor(df_bins$seqnames, levels = chrom_order)
df_bins <- df_bins[order(df_bins$seqnames, df_bins$start, df_bins$end), ]
df_bins$seqnames <- as.character(df_bins$seqnames)
row.names(df_bins) <- with(df_bins, paste(seqnames, start, end, sep = "_"))

list_bed_eigen <- list()

for (name in unlist(cellular_types_list)){
  print(name)
  
  bed <- read.table(paste0(bedgraph_eigen_dir, name, "_combined_compartments_eigen.bedgraph"), skip = 1)
  colnames(bed) <- c("seqnames", "start", "end", paste0("eigen_", name))
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  
  
  if (min(start(bed)) == 2) start(bed) <- start(bed) - 1
  
  list_bed_eigen[[name]] <- bed
}

names(list_bed_eigen) <- unlist(cellular_types_list)

df_bins_eigen <- as.data.frame(list_bed_eigen[['NPCwt_n1']])[c("seqnames", "start", "end")]

for (name in names(list_bed_eigen)){
  df_new_sample <- as.data.frame(list_bed_eigen[[name]])[c("seqnames", "start", "end", paste0("eigen_", name))]
  df_bins_eigen <- merge(df_bins_eigen, df_new_sample, all = TRUE, by = c("seqnames", "start", "end"))
}

chrom_order <- paste0("chr", c(1:19, "X"))
df_bins_eigen$seqnames <- factor(df_bins_eigen$seqnames, levels = chrom_order)
df_bins_eigen <- df_bins_eigen[order(df_bins_eigen$seqnames, df_bins_eigen$start, df_bins_eigen$end), ]
df_bins_eigen$seqnames <- as.character(df_bins_eigen$seqnames)
row.names(df_bins_eigen) <- with(df_bins_eigen, paste(seqnames, start, end, sep = "_"))


head(df_bins_eigen)
tail(df_bins_eigen)

## PCA of compartments
df_eigen_forPCA <- df_bins_eigen[, grepl("eigen", colnames(df_bins_eigen))]
data.pca <- prcomp(t(df_eigen_forPCA))

sample_to_group <- stack(lapply(cellular_types_list, function(x) x))
colnames(sample_to_group) <- c("sample", "group")

df_for_pca <- data.frame(
  sample = sub("eigen_", "", rownames(t(df_eigen_forPCA)))
)
df_for_pca <- merge(df_for_pca, sample_to_group, by = "sample", all.x = TRUE)
df_for_pca <- df_for_pca[match(sub("eigen_", "", rownames(t(df_eigen_forPCA))), df_for_pca$sample), ]

print(df_for_pca)


options(repr.plot.width=16, repr.plot.height=6)
#pdf("//storage-daredevil/sammyseq_nfcore/Analisi/MEFs/2025_11_25_compartments/compartments/consensus_maggioranza/pca_plots/pca_compartments.pdf",width = 8,height = 6)

fviz_pca_ind(data.pca#, label="none"
             #                 , habillage=df_for_pca$sample
             , addEllipses=FALSE, 
             geom.ind = c("text"),
             ellipse.level=0.95) + geom_point(aes(#shape = factor(df_for_pca$type), 
               color = factor(df_for_pca$sample),size=4,alpha = .5)) +
  #  scale_color_manual(values = c("#FDAC16","#57B5FB","#35AD78"))+
  theme_classic()  +
  theme(#text = element_text(family = "Arial") ,
    #                             axis.text=element_text(size=20),
    #                             plot.title = element_text(size=30),
    #                             legend.title = element_text(size=20),
    #                             legend.text = element_text(size=20),
    #                             axis.title=element_text(size=30),
    axis.line=element_blank(),
    #       axis.text.x=element_blank(),
    #       axis.text.y=element_blank(),
    #       axis.ticks=element_blank(),
    #                           axis.title.x=element_blank(),
    #                           axis.title.y=element_blank(),
    #                           legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()#,
    #                             panel.grid.major.x = element_line(color = "black",
    #                                           size = 0.1,
    #                                           linetype = 2)
  ) 
# dev.off()

## Consensus import
df_bins_forperc<-df_bins
df_bins_forperc$startend <-paste0(df_bins_forperc$seqnames,"_",df_bins_forperc$start,"_",df_bins_forperc$end)
rownames(df_bins_forperc)<-df_bins_forperc$startend
cons_path <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/compartments/consensus"

consensus_files <- list(
  MUT_NPC_consensus = "MUT_NPC_consensus_majority.bed",
  WT_NPC_consensus = "WT_NPC_consensus_majority.bed"
)

# Import each consensus file and add as column
for (col_name in names(consensus_files)) {
  file_path <- file.path(cons_path, consensus_files[[col_name]])
  
  bed <- read.table(file_path, 
                    header = FALSE, 
                    stringsAsFactors = FALSE,
                    skip = 1,  # Skip the "track name=" line
                    col.names = c("chr", "start", "end", "compartment", 
                                  "score", "strand", "thickStart", "thickEnd", "itemRgb"))
  
  df_bins_forperc[[col_name]] <- bed$compartment
  
  cat("Imported", col_name, "with", nrow(bed), "rows\n")
}


head(df_bins_forperc)

table(df_bins_forperc$MUT_NPC_consensus)

table(df_bins_forperc$WT_NPC_consensus)

spot_diff<- function(df_bins_forperc,consesustoeval,name1){
  
  x<- with(df_bins_forperc, ifelse(df_bins_forperc[[name1]] == "A" & eval(as.name(consesustoeval))  == "A" , 'A',
                                   ifelse(df_bins_forperc[[name1]] == "B" & eval(as.name(consesustoeval))   == "A" , 'A->B',
                                          
                                          ifelse(df_bins_forperc[[name1]] == "B" & eval(as.name(consesustoeval))   == "B", 'B',
                                                 
                                                 ifelse(df_bins_forperc[[name1]] == "A" & eval(as.name(consesustoeval))   == "B" , "B->A",
                                                        
                                                        "BOH" )))))
  return(x)
}
#Singoli exh10d against eff10d
MEFwt_c <- as.data.frame(table(df_bins_forperc$WT_NPC_consensus))
MEFwt_c$name <- rep("WT_NPC_consensus", nrow(MEFwt_c))

names(df_bins_forperc)

df_bins_forperc$NPCwt_n1vsMEFwt_c_df_toplot <- spot_diff(df_bins_forperc, "WT_NPC_consensus", "NPCwt_n1")
df_bins_forperc$NPCwt_n2vsMEFwt_c_df_toplot <- spot_diff(df_bins_forperc, "WT_NPC_consensus", "NPCwt_n2")
df_bins_forperc$NPCwt_IRE_R3vsMEFwt_c_df_toplot <- spot_diff(df_bins_forperc, "WT_NPC_consensus", "NPCwt_IRE_R3")
NPCwt_n1vsMEFwt_c_df <- as.data.frame(table(df_bins_forperc$NPCwt_n1vsMEFwt_c_df_toplot)); NPCwt_n1vsMEFwt_c_df$name <- "NPCwt_n1"
NPCwt_n2vsMEFwt_c_df <- as.data.frame(table(df_bins_forperc$NPCwt_n2vsMEFwt_c_df_toplot)); NPCwt_n2vsMEFwt_c_df$name <- "NPCwt_n2"
NPCwt_IRE_R3vsMEFwt_c_df <- as.data.frame(table(df_bins_forperc$NPCwt_IRE_R3vsMEFwt_c_df_toplot)); NPCwt_IRE_R3vsMEFwt_c_df$name <- "NPCwt_IRE_R3"

freq_preperc <- rbind(
  MEFwt_c,
  NPCwt_n1vsMEFwt_c_df,
  NPCwt_n2vsMEFwt_c_df,
  NPCwt_IRE_R3vsMEFwt_c_df
)
freq_preperc$name <- factor(
  freq_preperc$name,
  levels = c(
    "NPCwt_IRE_R3",
    "NPCwt_n2",
    "NPCwt_n1",
    "WT_NPC_consensus" 
  ),
  ordered = TRUE
)

freq_preperc$Var1 <- factor(
  freq_preperc$Var1,
  levels = c('B','A->B','B->A','A','equal','BOH'),
  ordered = TRUE
)

freq_preperc <- with(freq_preperc, freq_preperc[order(name),])
freq_preperc <- with(freq_preperc, freq_preperc[order(Var1),])

freq_preperc$perc <- (freq_preperc$Freq / nrow(df_bins_forperc)) * 100

freq_preperc <- freq_preperc[!freq_preperc$Var1 %in% c("equal","BOH"), ]

x_for_test <- freq_preperc[freq_preperc$Var1 != "B" & freq_preperc$Var1 != "A", ]
aggregate(perc ~ name, x_for_test, sum)

freq_preperc<- with(freq_preperc, freq_preperc[order(name),])
freq_preperc<- with(freq_preperc, freq_preperc[order(Var1),])
freq_preperc$perc <- (freq_preperc$Freq/nrow(df_bins_forperc))*100
freq_preperc<- freq_preperc[!freq_preperc$Var1== "equal",]
freq_preperc<- freq_preperc[!freq_preperc$Var1== "BOH",]

#pdf("//storage-daredevil/sammyseq_nfcore/Analisi/MEFs/2025_11_25_compartments/compartments/consensus_maggioranza/stacked_plots/single_replicates/MEFwt.pdf",height = 8,width = 10)
options(repr.plot.width=14, repr.plot.height=6)
ggplot(freq_preperc, aes(fill=Var1, y=name, x=perc), position = "fill") + 
  geom_bar(position="stack", stat="identity",width = 0.85 ,colour="transparent",) +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values=c(  "#E0AA58", "#F4DC96" ,"#C3DACF","#5A958F","red"),)+
  xlab("")+
  ylab(label = "% regions")+
  theme(axis.text=element_text(size=16,angle = 90),
        axis.title=element_text(size=14))+
  theme_classic()  +theme(
    axis.line=element_blank(),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()#,
  )
#dev.off()
names(df_bins_forperc)

#Singoli exh10d against eff10d
MEFACKO_c <- as.data.frame(table(df_bins_forperc$MUT_NPC_consensus))
MEFACKO_c$name <- rep("MUT_NPC_consensus", nrow(MEFACKO_c))

df_bins_forperc$NPCmutR1_IrevsMEFACKO_c_df_toplot <- spot_diff(df_bins_forperc, "MUT_NPC_consensus", "NPCmutR1_Ire")
df_bins_forperc$NPCmutR2_IrevsMEFACKO_c_df_toplot <- spot_diff(df_bins_forperc, "MUT_NPC_consensus", "NPCmutR2_Ire")
df_bins_forperc$NPCmut_IRE_R3vsMEFACKO_c_df_toplot <- spot_diff(df_bins_forperc, "MUT_NPC_consensus", "NPCmut_IRE_R3")
NPCmutR1_IrevsMEFACKO_c_df <- as.data.frame(table(df_bins_forperc$NPCmutR1_IrevsMEFACKO_c_df_toplot)); NPCmutR1_IrevsMEFACKO_c_df$name <- "NPCmutR1_Ire"
NPCmutR2_IrevsMEFACKO_c_df <- as.data.frame(table(df_bins_forperc$NPCmutR2_IrevsMEFACKO_c_df_toplot)); NPCmutR2_IrevsMEFACKO_c_df$name <- "NPCmutR2_Ire"
NPCmut_IRE_R3vsMEFACKO_c_df <- as.data.frame(table(df_bins_forperc$NPCmut_IRE_R3vsMEFACKO_c_df_toplot)); NPCmut_IRE_R3vsMEFACKO_c_df$name <- "NPCmut_IRE_R3"

freq_preperc <- rbind(
  MEFACKO_c,
  NPCmutR1_IrevsMEFACKO_c_df,
  NPCmutR2_IrevsMEFACKO_c_df,
  NPCmut_IRE_R3vsMEFACKO_c_df
)
freq_preperc$name <- factor(
  freq_preperc$name,
  levels = c(
    "NPCmutR1_Ire",
    "NPCmutR2_Ire",
    "NPCmut_IRE_R3",
    "MUT_NPC_consensus" 
  ),
  ordered = TRUE
)

freq_preperc$Var1 <- factor(
  freq_preperc$Var1,
  levels = c('B','A->B','B->A','A','equal','BOH'),
  ordered = TRUE
)

freq_preperc <- with(freq_preperc, freq_preperc[order(name),])
freq_preperc <- with(freq_preperc, freq_preperc[order(Var1),])

freq_preperc$perc <- (freq_preperc$Freq / nrow(df_bins_forperc)) * 100

freq_preperc <- freq_preperc[!freq_preperc$Var1 %in% c("equal","BOH"), ]

x_for_test <- freq_preperc[freq_preperc$Var1 != "B" & freq_preperc$Var1 != "A", ]
aggregate(perc ~ name, x_for_test, sum)

freq_preperc<- with(freq_preperc, freq_preperc[order(name),])
freq_preperc<- with(freq_preperc, freq_preperc[order(Var1),])
freq_preperc$perc <- (freq_preperc$Freq/nrow(df_bins_forperc))*100
freq_preperc<- freq_preperc[!freq_preperc$Var1== "equal",]
freq_preperc<- freq_preperc[!freq_preperc$Var1== "BOH",]

#pdf("//storage-daredevil/sammyseq_nfcore/Analisi/MEFs/2025_11_25_compartments/compartments/consensus_maggioranza/stacked_plots/single_replicates/MEFACKO.pdf",height = 8,width = 10)
options(repr.plot.width=14, repr.plot.height=6)
ggplot(freq_preperc, aes(fill=Var1, y=name, x=perc), position = "fill") + 
  geom_bar(position="stack", stat="identity",width = 0.85 ,colour="transparent",) +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values=c(  "#E0AA58", "#F4DC96" ,"#C3DACF","#5A958F","red"),)+
  xlab("")+
  ylab(label = "% regions")+
  theme(axis.text=element_text(size=16,angle = 90),
        axis.title=element_text(size=14))+
  theme_classic()  +theme(
    axis.line=element_blank(),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()#,
  )
# dev.off()

##singoli vs consenus contraria


MEFwt_c <- as.data.frame(table(df_bins_forperc$WT_NPC_consensus))
MEFwt_c$name <- rep("WT_NPC_consensus", nrow(MEFwt_c))

names(df_bins_forperc)

df_bins_forperc$NPCmutR1_IrevsMEFwt_c_df_toplot <- spot_diff(df_bins_forperc, "WT_NPC_consensus", "NPCmutR1_Ire")
df_bins_forperc$NPCmutR2_IrevsMEFwt_c_df_toplot <- spot_diff(df_bins_forperc, "WT_NPC_consensus", "NPCmutR2_Ire")
df_bins_forperc$NPCmut_IRE_R3vsMEFwt_c_df_toplot <- spot_diff(df_bins_forperc, "WT_NPC_consensus", "NPCmut_IRE_R3")
NPCmutR1_IrevsMEFwt_c_df <- as.data.frame(table(df_bins_forperc$NPCmutR1_IrevsMEFwt_c_df_toplot)); NPCmutR1_IrevsMEFwt_c_df$name <- "NPCmutR1_Ire"
NPCmutR2_IrevsMEFwt_c_df <- as.data.frame(table(df_bins_forperc$NPCmutR2_IrevsMEFwt_c_df_toplot)); NPCmutR2_IrevsMEFwt_c_df$name <- "NPCmutR2_Ire"
NPCmut_IRE_R3vsMEFwt_c_df <- as.data.frame(table(df_bins_forperc$NPCmut_IRE_R3vsMEFwt_c_df_toplot)); NPCmut_IRE_R3vsMEFwt_c_df$name <- "NPCmut_IRE_R3"


freq_preperc <- rbind(
  MEFwt_c,
  NPCmutR1_IrevsMEFwt_c_df,
  NPCmutR2_IrevsMEFwt_c_df,
  NPCmut_IRE_R3vsMEFwt_c_df
)
freq_preperc$name <- factor(
  freq_preperc$name,
  levels = c(
    "NPCmutR1_Ire",
    "NPCmutR2_Ire",
    "NPCmut_IRE_R3",
    "WT_NPC_consensus" 
  ),
  ordered = TRUE
)

freq_preperc$Var1 <- factor(
  freq_preperc$Var1,
  levels = c('B','A->B','B->A','A','equal','BOH'),
  ordered = TRUE
)

freq_preperc <- with(freq_preperc, freq_preperc[order(name),])
freq_preperc <- with(freq_preperc, freq_preperc[order(Var1),])

freq_preperc$perc <- (freq_preperc$Freq / nrow(df_bins_forperc)) * 100

freq_preperc <- freq_preperc[!freq_preperc$Var1 %in% c("equal","BOH"), ]

x_for_test <- freq_preperc[freq_preperc$Var1 != "B" & freq_preperc$Var1 != "A", ]
aggregate(perc ~ name, x_for_test, sum)

freq_preperc<- with(freq_preperc, freq_preperc[order(name),])
freq_preperc<- with(freq_preperc, freq_preperc[order(Var1),])
freq_preperc$perc <- (freq_preperc$Freq/nrow(df_bins_forperc))*100
freq_preperc<- freq_preperc[!freq_preperc$Var1== "equal",]
freq_preperc<- freq_preperc[!freq_preperc$Var1== "BOH",]

#pdf("//storage-daredevil/sammyseq_nfcore/Analisi/MEFs/2025_11_25_compartments/compartments/consensus_maggioranza/stacked_plots/single_replicates/MEFACKOvsMEFWT_single.pdf",height = 8,width = 10)
options(repr.plot.width=14, repr.plot.height=6)
ggplot(freq_preperc, aes(fill=Var1, y=name, x=perc), position = "fill") + 
  geom_bar(position="stack", stat="identity",width = 0.85 ,colour="transparent",) +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values=c(  "#E0AA58", "#F4DC96" ,"#C3DACF","#5A958F","red"),)+
  xlab("")+
  ylab(label = "% regions")+
  theme(axis.text=element_text(size=16,angle = 90),
        axis.title=element_text(size=14))+
  theme_classic()  +theme(
    axis.line=element_blank(),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()#,
  )
#dev.off()

## AtoB / BtoA switches
spot_diff<- function(df_bins_forperc,consesustoeval,name1){
  
  x<- with(df_bins_forperc, ifelse(df_bins_forperc[[name1]] == "A" & eval(as.name(consesustoeval))  == "A" , 'A',
                                   ifelse(df_bins_forperc[[name1]] == "B" & eval(as.name(consesustoeval))   == "A" , 'A->B',
                                          
                                          ifelse(df_bins_forperc[[name1]] == "B" & eval(as.name(consesustoeval))   == "B", 'B',
                                                 
                                                 ifelse(df_bins_forperc[[name1]] == "A" & eval(as.name(consesustoeval))   == "B" , "B->A",
                                                        
                                                        "BOH" )))))
  return(x)
}
spot_diff<- function(df_bins_forperc,consesustoeval,name1){
  
  x<- with(df_bins_forperc, ifelse(df_bins_forperc[[name1]] == "A" & eval(as.name(consesustoeval))  == "A" , 'A',
                                   ifelse(df_bins_forperc[[name1]] == "B" & eval(as.name(consesustoeval))   == "A" , 'A->B',
                                          
                                          ifelse(df_bins_forperc[[name1]] == "B" & eval(as.name(consesustoeval))   == "B", 'B',
                                                 
                                                 ifelse(df_bins_forperc[[name1]] == "A" & eval(as.name(consesustoeval))   == "B" , "B->A",
                                                        
                                                        "BOH" )))))
  return(x)
}

names(df_bins_forperc)
df_bins_forperc$NPCwt_n1_toplot<-spot_diff(df_bins_forperc,"WT_NPC_consensus","NPCwt_n1")
df_bins_forperc$NPCwt_n2_toplot<-spot_diff(df_bins_forperc,"WT_NPC_consensus","NPCwt_n2")
df_bins_forperc$NPCwt_IRE_R3_toplot<-spot_diff(df_bins_forperc,"WT_NPC_consensus","NPCwt_IRE_R3")

df_bins_forperc$NPCwt_n1CoV_toplot<-spot_diff(df_bins_forperc,"WT_NPC_consensus","NPCwt_n1")
df_bins_forperc$NPCwt_n2CoV_toplot<-spot_diff(df_bins_forperc,"WT_NPC_consensus","NPCwt_n2")
df_bins_forperc$NPCwt_IRE_R3CoV_toplot<-spot_diff(df_bins_forperc,"WT_NPC_consensus","NPCwt_IRE_R3")

MEFwt_c<-as.data.frame(table(df_bins_forperc$WT_NPC_consensus))
MEFwt_c$name<-rep("WT_NPC_consensus",times = nrow(MEFwt_c))
NPCwt_n1_df<-as.data.frame(table(df_bins_forperc$NPCwt_n1_toplot))
NPCwt_n1_df$name<-rep("NPCwt_n1",times = nrow(NPCwt_n1_df))
NPCwt_n2_df<-as.data.frame(table(df_bins_forperc$NPCwt_n2_toplot))
NPCwt_n2_df$name<-rep("NPCwt_n2",times = nrow(NPCwt_n2_df))
NPCwt_IRE_R3_df<-as.data.frame(table(df_bins_forperc$NPCwt_IRE_R3_toplot))
NPCwt_IRE_R3_df$name<-rep("NPCwt_IRE_R3",times = nrow(NPCwt_IRE_R3_df))

names(df_bins_forperc)

MEFACKO_c<-as.data.frame(table(df_bins_forperc$MUT_NPC_consensus))
MEFACKO_c$name<-rep("MUT_NPC_consensus",times = nrow(MEFACKO_c))

NPCwt_n1CoV_df<-as.data.frame(table(df_bins_forperc$NPCwt_n1CoV_toplot))
NPCwt_n1CoV_df$name<-rep("NPCmutR1_Ire",times = nrow(NPCwt_n1CoV_df))
NPCwt_n2CoV_df<-as.data.frame(table(df_bins_forperc$NPCwt_n2CoV_toplot))
NPCwt_n2CoV_df$name<-rep("NPCmutR2_Ire",times = nrow(NPCwt_n2CoV_df))
NPCwt_IRE_R3CoV_df<-as.data.frame(table(df_bins_forperc$NPCwt_IRE_R3CoV_toplot))
NPCwt_IRE_R3CoV_df$name<-rep("NPCmut_IRE_R3",times = nrow(NPCwt_IRE_R3CoV_df))
freq_preperc<- rbind(MEFwt_c,NPCwt_n1_df,NPCwt_n2_df,NPCwt_IRE_R3_df,MEFACKO_c,NPCwt_n1CoV_df,NPCwt_n2CoV_df,NPCwt_IRE_R3CoV_df)

freq_preperc$name<-factor(freq_preperc$name,levels = c("NPCmut_IRE_R3","NPCmutR2_Ire","NPCmutR1_Ire","MUT_NPC_consensus","NPCwt_IRE_R3","NPCwt_n2","NPCwt_n1","WT_NPC_consensus"),ordered = TRUE)
freq_preperc$Var1<-factor(freq_preperc$Var1,levels = c('B',"B->A",'A->B','A',"equal","BOH"),ordered = TRUE)

freq_preperc <- with(freq_preperc, freq_preperc[order(name),])
freq_preperc <- with(freq_preperc, freq_preperc[order(Var1),])

freq_preperc$perc <- (freq_preperc$Freq / nrow(df_bins_forperc)) * 100

freq_preperc <- freq_preperc[!freq_preperc$Var1 %in% c("equal","BOH"), ]

x_for_test <- freq_preperc[freq_preperc$Var1 != "B" & freq_preperc$Var1 != "A", ]
aggregate(perc ~ name, x_for_test, sum)

freq_preperc<- with(freq_preperc, freq_preperc[order(name),])
freq_preperc<- with(freq_preperc, freq_preperc[order(Var1),])
freq_preperc$perc <- (freq_preperc$Freq/nrow(df_bins_forperc))*100
freq_preperc<- freq_preperc[!freq_preperc$Var1== "equal",]
freq_preperc<- freq_preperc[!freq_preperc$Var1== "BOH",]

#pdf("//storage-daredevil/sammyseq_nfcore/Analisi/MEFs/2025_11_25_compartments/compartments/consensus_maggioranza/stacked_plots/MEFACKOVSMEFwt_model",height = 8,width = 10)
options(repr.plot.width=14, repr.plot.height=6)
ggplot(freq_preperc, aes(fill=Var1, y=name, x=perc), position = "fill") + 
  geom_bar(position="stack", stat="identity",width = 0.85 ,colour="transparent",) +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  #position="fill" "stack"
  #     scale_fill_viridis(discrete = T) +
  scale_fill_manual(values=c(  "#E0AA58", "#F4DC96" ,"#C3DACF","#5A958F","red"),)+
  xlab("")+
  ylab(label = "% regions")+
  theme(axis.text=element_text(size=16,angle = 90),
        axis.title=element_text(size=14))+
  
  #   facet_grid( scales = "free_y", space = "free_x") + #ylim(c(0,8))+
  # ggtitle("cpcdde")+
  theme_classic()  +theme(#text = element_text(family = "Arial") ,
    #                             axis.text=element_text(size=12),
    #                             plot.title = element_text(size=12),
    legend.title = element_blank(),
    #                             legend.text = element_text(size=12),
    #                             axis.title=element_text(size=12),
    axis.line=element_blank(),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    #       axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    #                           legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    #                             panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()#,
    #                             panel.grid.major.x = element_line(color = "black",
    #                                           size = 0.1,
    #                                           linetype = 2)
  )
#dev.off()

## Comparison between consensus groups
df_bins_foralluvial<-df_bins_forperc

df_bins_foralluvial$MUTNPCvsWTNPC_toplot<-spot_diff(df_bins_foralluvial,"WT_NPC_consensus","MUT_NPC_consensus")

table(df_bins_foralluvial$MUTNPCvsWTNPC_toplot)

df_forrich<- df_bins_foralluvial[,c("WT_NPC_consensus","MUTNPCvsWTNPC_toplot")]

wtcc<-as.data.frame(table(df_forrich$WT_NPC_consensus))
wtcc$name<-rep("WT_NPC_consensus",times = nrow(wtcc))

b11<-as.data.frame(table(df_forrich$MUTNPCvsWTNPC_toplot))
b11$name<-rep("MUTNPCvsWTNPC_toplot",times = nrow(b11))

freq_preperccc<- rbind(wtcc,b11)

freq_preperccc$name<-factor(freq_preperccc$name,levels = c("MUTNPCvsWTNPC_toplot","WT_NPC_consensus"),ordered = TRUE)
freq_preperccc$Var1<-factor(freq_preperccc$Var1,levels = c('B',"B->A",'A->B','A',"BOH","equal"),ordered = TRUE)

freq_preperccc<- with(freq_preperccc, freq_preperccc[order(name),])
freq_preperccc<- with(freq_preperccc, freq_preperccc[order(Var1),])
freq_preperccc$perc <- (freq_preperccc$Freq/nrow(df_bins_forperc))*100


freq_preperccc$Var1<-factor(freq_preperccc$Var1,levels = c('B','A->B',"B->A",'A',"BOH","equal"),ordered = TRUE)

freq_preperccc<- with(freq_preperccc, freq_preperccc[order(name),])
freq_preperccc<- with(freq_preperccc, freq_preperccc[order(Var1),])
freq_preperccc$perc <- (freq_preperccc$Freq/nrow(df_bins_forperc))*100

freq_preperccc<- freq_preperccc[!freq_preperccc$Var1== "BOH",]
freq_preperccc<- freq_preperccc[!freq_preperccc$Var1== "equal",]
# aggregate(width ~ seqnames  + bintype,  df_forrich,sum)
#pdf("//storage-daredevil/sammyseq_nfcore/Analisi/MEFs/2025_11_25_compartments/compartments/consensus_maggioranza/stacked_plots/MEFACKOvsMEFwt_model",height = 8,width = 10)
options(repr.plot.width=16, repr.plot.height=6)
ggplot(freq_preperccc, aes(fill=Var1, y=name, x=perc), position = "fill") + 
  geom_bar(position="stack", stat="identity",width = 0.85 ,colour="transparent",) +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  #position="fill" "stack"
  #     scale_fill_viridis(discrete = T) +
  scale_fill_manual(values=c(  "#E0AA58", "#F4DC96" ,"#C3DACF","#5A958F"),)+
  xlab("")+
  ylab(label = "% regions")+
  theme(axis.text=element_text(size=16,angle = 90),
        axis.title=element_text(size=14))+
  
  #   facet_grid( scales = "free_y", space = "free_x") + #ylim(c(0,8))+
  # ggtitle("cpcdde")+
  theme_classic()  +theme(#text = element_text(family = "Arial") ,
    #                             axis.text=element_text(size=12),
    #                             plot.title = element_text(size=12),
    legend.title = element_blank(),
    #                             legend.text = element_text(size=12),
    #                             axis.title=element_text(size=12),
    axis.line=element_blank(),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    #       axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    #                           legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    #                             panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()#,
    #                             panel.grid.major.x = element_line(color = "black",
    #                                           size = 0.1,
    #                                           linetype = 2)
  )
#dev.off()


txdb_name<- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/gencode.vM25.basic.annotation.gtf"
txdb <- GenomicFeatures::makeTxDbFromGFF(txdb_name)
genes <- GenomicFeatures::genes(txdb)

summarizeProteinCodingGenes <- function(txdb){
  stopifnot(is(txdb, "TxDb"))
  protein_coding_tx <- names(GenomicFeatures::cdsBy(txdb, use.names=TRUE))
  all_tx <- mcols(GenomicFeatures::transcripts(txdb, columns=c("gene_id", "tx_name")))
  all_tx$gene_id <- as.character(all_tx$gene_id)
  all_tx$is_coding <- all_tx$tx_name %in% protein_coding_tx
  tmp <- splitAsList(all_tx$is_coding, all_tx$gene_id)
  gene <- names(tmp)
  n_tx <- unname(sum(tmp, na.rm = TRUE))
  n_coding <- unname(sum(tmp))
  n_non_coding <- n_tx - n_coding
  data.frame(gene, n_tx, n_coding, n_non_coding, stringsAsFactors=FALSE)
}

geneid_codingdf  <- summarizeProteinCodingGenes(txdb)
genes_gr         <- genes[genes$gene_id %in% geneid_codingdf[geneid_codingdf$n_coding > 0, ]$gene]


genes_gr$gene_id<-gsub("\\..*","",genes_gr$gene_id)

# BTOA
MUTNPCvsWTNPC_toplot_ALL_BTOA<-df_bins_foralluvial[df_bins_foralluvial$MUTNPCvsWTNPC_toplot=="B->A",]
MUTNPCvsWTNPC_toplot_ALL_BTOA_gr<- makeGRangesFromDataFrame(MUTNPCvsWTNPC_toplot_ALL_BTOA,keep.extra.columns = TRUE)
write.table(genes_gr[findOverlaps(MUTNPCvsWTNPC_toplot_ALL_BTOA_gr,promoters(genes_gr,upstream = 2500,downstream = 500))@to]$gene_id,
            paste0("/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/compartments/consensus/gene_id_","MUTNPCvsWTNPC_toplot_ALL_BTOA_gr",".txt"),
            quote = FALSE,
            sep="\t",
            row.names=FALSE,
            col.names = TRUE
)

# ATOB
MUTNPCvsWTNPC_toplot_ALL_ATOB<-df_bins_foralluvial[df_bins_foralluvial$MUTNPCvsWTNPC_toplot=="A->B",]
MUTNPCvsWTNPC_toplot_ALL_ATOB_gr<- makeGRangesFromDataFrame(MUTNPCvsWTNPC_toplot_ALL_ATOB,keep.extra.columns = TRUE)
write.table(genes_gr[findOverlaps(MUTNPCvsWTNPC_toplot_ALL_ATOB_gr,promoters(genes_gr,upstream = 2500,downstream = 500))@to]$gene_id,
            paste0("/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/compartments/consensus/gene_id_","MUTNPCvsWTNPC_toplot_ALL_ATOB_gr",".txt"),
            quote = FALSE,
            sep="\t",
            row.names=FALSE,
            col.names = TRUE
)

library(GenomicRanges)
library(GenomicFeatures)

# load your bins (A/B for WT or MUT)
bins <- makeGRangesFromDataFrame(df_bins_forperc, keep.extra.columns = TRUE)

# load gene annotation
txdb <- makeTxDbFromGFF("gencode.vM25.annotation.gtf")  # <- adjust accordingly
genes <- genes(txdb)

# find genes overlapping any compartment bin
hits <- findOverlaps(genes, bins)
background_genes <- unique(genes$gene_id[queryHits(hits)])
# Remove the version number after the dot
background_genes_clean <- unique(sub("\\..*$", "", background_genes))

writeLines(background_genes_clean, "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/compartments_analysis/compartments/consensus/background_genes.txt")
