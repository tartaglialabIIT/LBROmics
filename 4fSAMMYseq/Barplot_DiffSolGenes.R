library(tidyverse)

plot_gene_counts <- function(folder, pattern, pdf_file) {
  
  # --- List all files matching pattern ---
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if(length(files) == 0) stop("No files found with this pattern!")
  
  # --- Function to parse each file ---
  parse_filename <- function(fname) {
    base <- basename(fname)
    parts <- strsplit(base, "_")[[1]]
    
    # Detect condition
    condition <- paste(parts[2:min(6, length(parts)-3)], collapse = "_")
    
    # Extract stage and direction (last 2 parts before file extension)
    stage_dir <- parts[length(parts)-2]
    direction <- parts[length(parts)-1]
    category <- paste0(stage_dir, "_", direction)
    
    genes <- readLines(fname)
    genes <- genes[genes != ""]
    
    tibble(
      file = base,
      condition = condition,
      category = category,
      n_genes = length(genes)
    )
  }
  
  # --- Apply parsing ---
  df <- map_dfr(files, parse_filename)
  
  # --- Clean condition names ---
  df$condition <- case_when(
    grepl("WT_NPC_vs_WT_ESC", df$condition) ~ "WT_NPC_vs_WT_ESC",
    grepl("MUT_NPC_vs_WT_NPC", df$condition) ~ "MUT_NPC_vs_WT_NPC",
    grepl("MUT_ESC_vs_WT_ESC", df$condition) ~ "MUT_ESC_vs_WT_ESC",
    TRUE ~ df$condition
  )
  df <- df %>% drop_na(condition)
  
  # --- Define stages and categories in order ---
  stages <- strsplit(pattern, "vs")[[1]]
  categories <- c(
    paste0(stages[1], "_up"),
    paste0(stages[1], "_down"),
    paste0(stages[2], "_up"),
    paste0(stages[2], "_down")
  )
  
  # --- Define colors for these categories ---
  nice_cols <- setNames(c("red","orange","lightblue","darkblue"), categories)
  
  # --- Factor ordering ---
  df$category <- factor(df$category, levels = categories)
  df$condition <- factor(df$condition, levels = unique(df$condition))
  
  # --- Custom x-axis labels ---
  custom_labels <- c(
    "MUT_ESC_vs_WT_ESC" = "Lbr NT-KO vs WT (ESC)",
    "MUT_NPC_vs_WT_NPC" = "Lbr NT-KO vs WT (NPC)",
    "WT_NPC_vs_WT_ESC"  = "WT NPC vs WT ESC"
  )
  
  # --- Create plot ---
  p <- ggplot(df, aes(x = condition, y = n_genes, fill = category)) +
    geom_col(width = 0.35, colour = "black") +
    geom_text(
      data = df %>% filter(condition %in% c("MUT_NPC_vs_WT_NPC", "WT_NPC_vs_WT_ESC")),
      aes(label = n_genes),
      position = position_stack(vjust = 0.5),
      size = 6, colour = "white", fontface = "bold"
    ) +
    scale_fill_manual(values = nice_cols) +
    scale_x_discrete(labels = custom_labels) +
    labs(x = NULL, y = "Number of genes", fill = NULL) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 20),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # --- Save as PDF ---
  ggsave(pdf_file, p, width = 8, height = 10)  # reduced width, increased height
  
  return(p)
}
# --- Example usage ---
folder <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/differential_solubility/genes/"
wd <- "/Users/jfiorentino/Desktop/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT/Cerase_single_cell/SAMMY-seq/differential_solubility/"
setwd(wd)

plot_gene_counts(folder = folder, pattern = "S2SvsS3", pdf_file = "../FINAL_PLOTS_NOV2025/diffsoluble_gene_counts_S2SvsS3.pdf")