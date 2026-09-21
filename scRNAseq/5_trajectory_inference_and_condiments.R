#!/usr/bin/env Rscript
# Note that this analysis runs in R4.2.2
library(Seurat)
library(slingshot)
library(tradeSeq)
library(UpSetR)
library(scales)
library(viridis)
library(condiments)

# For data manipulation
library(dplyr)
library(tidyr)

# For visualization
library(ggplot2)
library(RColorBrewer)
library(viridis)

theme_set(theme_classic())

# Set LBROMICS_SCRNA_ROOT to the Seurat analysis directory.
wd <- Sys.getenv("LBROMICS_SCRNA_ROOT", unset = "")
if (!nzchar(wd)) wd <- getwd()
setwd(wd)


# Load the clustered rds object with the diffusion map computed with destiny
data.combined <- readRDS("./files/sub_data_combined_and_clustered_dm.rds")
data.combined$cond <- factor(data.combined$cond, levels = c('WT', 'mutant'))
data.combined@meta.data$integrated_snn_res.0.2 <- factor(x = data.combined@meta.data$integrated_snn_res.0.2
                                                         ,levels=c(0,2,3,1))
data.combined@meta.data$integrated_snn_res.0.3 <- factor(x = data.combined@meta.data$integrated_snn_res.0.3
                                                         ,levels=c(1,2,3,4,0,5))

dimred <- data.combined@reductions$dm@cell.embeddings
clustering <- data.combined$integrated_snn_res.0.3
sds<- slingshot(dimred,clustering,start.clus="1")

curves <- slingCurves(sds, as.df = TRUE)

df <- data.frame(dimred, "Cluster" = as.character(clustering))
df$Cluster <- factor(df$Cluster,levels = c("1","2","3","4","0","5"))

pal <- hue_pal()(7)
pal <- c(pal[1],pal[2],pal[3],pal[5],pal[6],pal[7])

pdf("./plots_JULY_2023/trajectory_inference/slingshot_curves.pdf",
    width = 7, height=4)
p <- ggplot(df, aes(x = DC_1, y = DC_2)) +
  geom_point(aes(fill = Cluster), col = "grey70", shape = 21) + 
  scale_fill_manual(values=pal) + 
  theme_classic()
p + geom_path(data = curves %>% arrange(Order),
              aes(group = Lineage, col = as.character(Lineage)), size = 1.5) 
dev.off()
library(rgl)
library(rglwidget)
options(rgl.printRglwidget = TRUE)

pdf("./plots_JULY_2023/trajectory_inference/3d_trajectory.pdf")
plot3d(dimred, col = 'grey50', aspect = 'iso')
plot3d.SlingshotDataSet(SlingshotDataSet(sds),lwd=3, type = 'curves',add=T)
dev.off()
# Load the lineages from slingshot
load("./files_JULY_2023/curves.RData")

set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        #end.clus = c("11","7","10","9","5"), #define how many branches/lineages to consider
                        start.clus = "1") #define where to start the trajectories

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
save(curves, file = "./files_JULY_2023/curves.RData")

df <- data.frame(dimred)
df$cond <- data.combined$cond
df$cl <- data.combined$integrated_snn_res.0.3

top_res <- topologyTest(sds = curves, conditions = df$cond)
knitr::kable(top_res)

dir.create("./plots_JULY_2023/condiments")
write.csv(top_res,"./plots_JULY_2023/condiments/topologyTest.csv")

# The test says we can infer a common trajectory between mutant and WT

# Compute an imbalance score
scores <- imbalance_score(Object = df %>% select(DC_1, DC_2) %>% as.matrix(),
                          conditions = df$cond)
df$scores <- scores$scores
df$scaled_scores <- scores$scaled_scores

pdf(file = "./plots_JULY_2023/condiments/imbalance_score.pdf",width=8, height=6)
ggplot(df, aes(x = DC_1, y = DC_2, col = scores)) +
  geom_point() +
  scale_color_viridis_c(option = "C")
dev.off()

pdf(file = "./plots_JULY_2023/condiments/scaled_imbalance_score.pdf",width=8, height=6)
ggplot(df, aes(x = DC_1, y = DC_2, col = scaled_scores)) +
  geom_point() +
  scale_color_viridis_c(option = "C")
dev.off()

pt <- data.frame(slingPseudotime(curves))

psts <- slingPseudotime(curves) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = df$cond) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

pdf(file = "./plots_JULY_2023/condiments/diff_progression.pdf",width=10, height=5)
ggplot(psts, aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values=c('#F8766D','#00BFC4'))+
#  scale_fill_brewer(type = "qual") +
  facet_wrap(~lineages) +
  theme(legend.position = "bottom")
dev.off()

# Testing for differential progression 
prog_res <- progressionTest(curves, conditions = df$cond, global = TRUE, lineages = TRUE)
knitr::kable(prog_res)
write.csv(prog_res,"./plots_JULY_2023/condiments/diffProg.csv")

# Differential differentiation
# Lineage 1
df$weight_1 <- slingCurveWeights(curves, as.probs = TRUE)[, 1]
pdf("./plots_JULY_2023/condiments/diff_diff_lineage1.pdf",width=6, height=5)
ggplot(df, aes(x = weight_1, fill = cond)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values=c('#F8766D','#00BFC4'))+
#  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for the first lineage")
dev.off()

# Lineage 2
df$weight_2 <- slingCurveWeights(curves, as.probs = TRUE)[, 2]
pdf("./plots_JULY_2023/condiments/diff_diff_lineage2.pdf",width=6, height=5)
ggplot(df, aes(x = weight_2, fill = cond)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values=c('#F8766D','#00BFC4'))+
#  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for the second lineage")
dev.off()

# Differential differentiation test
set.seed(12)
dif_res <- differentiationTest(curves, conditions = df$cond, global = TRUE, pairwise = TRUE)
knitr::kable(dif_res)
write.csv(dif_res,"./plots_JULY_2023/condiments/diffDiff.csv")

counts <- as.matrix(data.combined@assays$RNA@counts)

# Choose number of knots for tradeSeq analysis
icMat <- evaluateK(counts = counts, sds = curves, k = 3:10, 
                   nGenes = 200, verbose = T,plot=TRUE)

