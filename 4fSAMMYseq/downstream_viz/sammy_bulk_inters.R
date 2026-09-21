library(tidyverse)
library(data.table)
library(ggplot2)
library(readr)

wd <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/"
setwd(wd)

library(tidyverse)
library(data.table)
library(ggplot2)
library(readr)

plot_sammy_bulk_overlap <- function(
    comparison = "S2SvsS3",
    group = "MUT_ESCvsWT_ESC",
    bulk_path,
    sammy_rds_path,
    save_pdf = NULL    # optional PDF output
) {
  
  ### --- 1. Load and classify bulk RNA-seq DEGs --- ###
  bulk <- fread(bulk_path) %>%
    mutate(direction = case_when(
      padj < 0.05 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & log2FoldChange < -0.5 ~ "DOWN",
      TRUE ~ "NS"
    ))
  
  bulk_up   <- bulk$gene_name[bulk$direction == "UP"]
  bulk_down <- bulk$gene_name[bulk$direction == "DOWN"]
  print(c(length(bulk_up),length(bulk_down)))
  all_bulk_genes <- bulk$gene_name
  
  print(length(all_bulk_genes))
  
  
  ### --- 2. Load SAMMY results --- ###
  sammy <- readRDS(sammy_rds_path)
  
  prefix <- gsub("vs", "_vs_", group)
  
  ### --- 2b. Build SAMMY lists dynamically based on comparison --- ###
  
  # Example comparison: "S2SvsS3"
  sides <- unlist(strsplit(comparison, "vs"))
  side1 <- sides[1]
  side2 <- sides[2]
  
  # category names (for plotting)
  categories <- c(
    paste0(side1, "_up"),
    paste0(side1, "_down"),
    paste0(side2, "_up"),
    paste0(side2, "_down")
  )
  
  # build internal SAMMY names based on group (e.g. MUT_ESC_vs_WT_ESC)
  sammy_internal_names <- c(
    paste0(prefix, "_", side1, "_up"),
    paste0(prefix, "_", side1, "_down"),
    paste0(prefix, "_", side2, "_up"),
    paste0(prefix, "_", side2, "_down")
  )
  
  # build the actual list
  sammy_lists <- purrr::map(sammy_internal_names, ~ sammy$genes[[.x]])
  
  # attach category names
  names(sammy_lists) <- categories
  
  ### --- 3. Compute universe = union of all genes seen anywhere --- ###
  all_sammy_genes <- unique(unlist(sammy_lists))
  universe <- union(all_bulk_genes, all_sammy_genes)
  universe_size <- length(universe)
  
  print(universe_size)
  
  ### --- 4. Compute intersections + Fisher tests --- ###
  results <- map_df(names(sammy_lists), function(cat) {
    
    genes <- sammy_lists[[cat]]
    
    # intersections
    int_up   <- length(intersect(genes, bulk_up))
    int_down <- length(intersect(genes, bulk_down))
    
    # fisher test using proper universe
    fisher_p <- function(intersect_n, bulk_vec) {
      
      a <- intersect_n
      b <- length(genes) - intersect_n
      c <- length(bulk_vec) - intersect_n
      d <- universe_size - a - b - c
      
      # ensure no negative cells (can happen rarely if vectors overlap oddly)
      d <- max(d, 0)
      
      mat <- matrix(c(a,b,c,d), nrow = 2)
      fisher.test(mat)$p.value
    }
    
    tibble(
      category = cat,
      sammy_size = length(genes),
      bulk_type = c("Bulk_UP", "Bulk_DOWN"),
      intersection = c(int_up, int_down),
      p = c(fisher_p(int_up, bulk_up), fisher_p(int_down, bulk_down))
    )
  })
  
  
  ### --- 5. Add significance asterisks --- ###
  results <- results %>%
    mutate(star = ifelse(p < 0.05, "*", "")) %>%
    mutate(category_label = paste0(category, " (N=", sammy_size, ")"))
  
  print(results)
  ### --- 6. Plot --- ###
  p <- ggplot(results, aes(x = category_label, y = intersection, fill = bulk_type)) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_text(aes(label = star),
              position = position_dodge(width = 0.8),
              vjust = -0.3,
              size = 6,
              fontface = "bold") +
    theme_bw(base_size = 14) +
    labs(
      y = "Intersection size",
      x = "SAMMY category",
      fill = "Bulk DEG type",
      title = paste("Overlap between SAMMY", comparison, "and bulk", group)
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ### --- 7. Save output to PDF if requested --- ###
  if (!is.null(save_pdf)) {
    ggsave(save_pdf, p, width = 8.5, height = 5.5)
  }
  
  return(p)
}

p <- plot_sammy_bulk_overlap(
  comparison = "S2SvsS3",
  group = "MUT_ESCvsWT_ESC",
  bulk_path = "./DESEQ2/DiffExp/deseq2_results_ESC.csv",
  sammy_rds_path = "./differential_solubility/rdata/S2SvsS3_MUT_ESCvsWT_ESC_analysis.rds",
  save_pdf = "overlap_ESC_S2SvsS3.pdf"
)

p   # prints the plot

p <- plot_sammy_bulk_overlap(
  comparison = "S2SvsS3",
  group = "MUT_NPCvsWT_NPC",
  bulk_path = "./DESEQ2/DiffExp/deseq2_results_NPC.csv",
  sammy_rds_path = "./differential_solubility/rdata/S2SvsS3_MUT_NPCvsWT_NPC_analysis.rds",
  save_pdf = "overlap_NPC_S2SvsS3.pdf"
)

p   # prints the plot

p <- plot_sammy_bulk_overlap(
  comparison = "S2LvsS3",
  group = "MUT_NPCvsWT_NPC",
  bulk_path = "./DESEQ2/DiffExp/deseq2_results_NPC.csv",
  sammy_rds_path = "./differential_solubility/rdata/S2LvsS3_MUT_NPCvsWT_NPC_analysis.rds",
  save_pdf = "overlap_NPC_S2LvsS3.pdf"
)

p   # prints the plot

p <- plot_sammy_bulk_overlap(
  comparison = "S2SvsS4",
  group = "MUT_NPCvsWT_NPC",
  bulk_path = "./DESEQ2/DiffExp/deseq2_results_NPC.csv",
  sammy_rds_path = "./differential_solubility/rdata/S2SvsS4_MUT_NPCvsWT_NPC_analysis.rds",
  save_pdf = "overlap_NPC_S2SvsS4.pdf"
)

p   # prints the plot

p <- plot_sammy_bulk_overlap(
  comparison = "S4vsS3",
  group = "MUT_NPCvsWT_NPC",
  bulk_path = "./DESEQ2/DiffExp/deseq2_results_NPC.csv",
  sammy_rds_path = "./differential_solubility/rdata/S4vsS3_MUT_NPCvsWT_NPC_analysis.rds",
  save_pdf = "overlap_NPC_S4vsS3.pdf"
)

p   # prints the plot
