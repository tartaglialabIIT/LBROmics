library(rlist)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(grid)
library(scales)
library(gprofiler2)
library(rrvgo)
library(htmlwidgets)

suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

add.col<-function(df, new.col){
  n.row<-dim(df)[1]
  length(new.col)<-n.row
  cbind(df, new.col)
}

FindSpecificMarkers <- function(seuratObject,idents.label,cl.label){
  DefaultAssay(seuratObject) <- "RNA"
  mito.genes <- grep(pattern = "^mt-", x = rownames(seuratObject), value = TRUE)
  #rps.genes <- grep(pattern = "^Rps", x = rownames(data.combined), value = TRUE)
  #rpl.genes <- grep(pattern = "^Rpl", x = rownames(data.combined), value = TRUE)
  
  genes.use <- rownames(seuratObject)[! rownames(seuratObject) %in% mito.genes]
  
  Idents(seuratObject) <- idents.label
  clusters<- data.frame(table(Idents(seuratObject)))[['Var1']]
  k.PC<-0
  for (j in 0:(length(clusters)-1)){
    if(cl.label!=j){
      tmp.unique.markers <- FindMarkers(seuratObject,
                                        ident.1 = cl.label,
                                        ident.2 = j,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25,
                                        genes.use = genes.use)
      #tmp.unique.markers <- tmp.unique.markers[(tmp.unique.markers$avg_logFC) > 1,]
      tmp.unique.markers <- tmp.unique.markers[(tmp.unique.markers$p_val_adj) < 0.05,]
      #print(length(rownames(tmp.unique.markers)))
      if(k.PC==0){
        unique.markers.df <- data.frame(rownames(tmp.unique.markers)
                                        , rep(1, length(rownames(tmp.unique.markers))),-log10(tmp.unique.markers$p_val_adj))
        colnames(unique.markers.df) <- c('Gene','Succ_PC','logp')
      }else{
        if(length(rownames(tmp.unique.markers))>0){
          for(k in 1:length(rownames(tmp.unique.markers))){
            if(rownames(tmp.unique.markers)[k] %in% unique.markers.df$Gene){
              idx <- match(rownames(tmp.unique.markers)[k],unique.markers.df$Gene)
              unique.markers.df$Succ_PC[idx] <- unique.markers.df$Succ_PC[idx]+1
              unique.markers.df$logp[idx] <- unique.markers.df$logp[idx]-log10(tmp.unique.markers$p_val_adj[k])
            }
          } 
        }
        
      } 
      #print(unique.markers.df)
      k.PC <- k.PC + 1
    }
  }
  unique.markers.df <- unique.markers.df[unique.markers.df$Succ_PC>=length(clusters)-2,]
  unique.markers.df <- unique.markers.df[order(-unique.markers.df$logp),]
  unique.markers.df
}

FindDECond <- function(seuratObject,idents.label){
  DefaultAssay(seuratObject) <- "RNA"
  mito.genes <- grep(pattern = "^mt-", x = rownames(seuratObject), value = TRUE)
  #rps.genes <- grep(pattern = "^Rps", x = rownames(data.combined), value = TRUE)
  #rpl.genes <- grep(pattern = "^Rpl", x = rownames(data.combined), value = TRUE)
  
  genes.use <- rownames(seuratObject)[! rownames(seuratObject) %in% mito.genes]
  
  Idents(seuratObject) <- idents.label
  clusters<- data.frame(table(Idents(seuratObject)))[['Var1']]
  seuratObject$celltype.condition <- paste(Idents(seuratObject), seuratObject$cond, sep="_")
  seuratObject$celltype <- Idents(seuratObject)
  Idents(seuratObject) <- "celltype.condition"
  my.list.up <- list()
  my.list.down <- list()
  for (i in 0:(length(clusters)-1)){
    ident1 <- paste0(i,"_WT")
    ident2 <- paste0(i,"_mutant")
    condition.diffgenes <- FindMarkers(seuratObject, ident.1 = ident1, ident.2=ident2,only.pos = FALSE, min.pct=0.25, logfc.threshold=0.25,genes.use = genes.use)
    condition.diffgenes <- condition.diffgenes[(condition.diffgenes$p_val_adj) < 0.05,]
    my.list.up[[i+1]] <- condition.diffgenes[(condition.diffgenes$avg_log2FC) > 0,]
    my.list.down[[i+1]] <- condition.diffgenes[(condition.diffgenes$avg_log2FC) < 0,]
  }
  list("up" = my.list.up, "down" = my.list.down)
}

custom.GO <- function(my.query,i,dir.id,my.bg,my.terms=FALSE){
  # Run GO enrichment analysis with gprofiler
  gprofiler_results <- gost(query = my.query, 
                            organism = "mmusculus", ordered_query = FALSE, 
                            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                            measure_underrepresentation = FALSE, evcodes = TRUE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "annotated", custom_bg = my.bg, 
                            numeric_ns = "", sources = NULL, as_short_link = FALSE)
  # Plot results
  p <- gostplot(gprofiler_results, capped = TRUE, interactive = T)
  
  my.Dir <- dir.id
  dir.create(file.path("./plots_JULY_2023/", my.Dir), showWarnings = FALSE)
  full.dir <- paste("./plots_JULY_2023",my.Dir,sep = '/')
  file.name <- paste("gprof",i,sep = "_")
  file.name <- paste(file.name,"html",sep=".")
  # Save the interactive plot to a html file
  saveWidget(p, file=paste(full.dir,file.name,sep = "/"))
  mydf <- as.data.frame(gprofiler_results$result)
  file.name <- paste("gprof_table",i,sep = "_")
  file.name <- paste(file.name,"csv",sep=".")
  file.name <- paste(full.dir,file.name,sep = "/")
  mydf <- apply(mydf,2,as.character)
  
  write.csv(x=mydf,file=file.name)
  
  if(my.terms!=FALSE){
    
    file.name <- paste("gprof",i,sep = "_")
    file.name <- paste(file.name,"pdf",sep=".")
    pdf(file = paste(full.dir,file.name,sep = "/"), width = 8,height = 6)
    # Plot for highlighting specific terms (should be run after gostplot 
    # (with interactive=FALSE))
    p <- gostplot(gprofiler_results, capped = TRUE, interactive = F)
    pp <- publish_gostplot(p, highlight_terms = my.terms, 
                             width = NA, height = NA, filename = NULL )
    pp
    dev.off()
  }
  
  
  # Extract the significant GO terms to be passed to REVIGO
  ## Subset and reorder gProfiler results to only include columns of interest
  gprofiler_results_reordered <- gprofiler_results$result[, c("term_id", "effective_domain_size", "term_name", "p_value", "intersection_size", "term_size", "intersection")]
  
  ## Order the results by p-adjusted value
  gprofiler_results_reordered <- gprofiler_results_reordered[order(gprofiler_results_reordered$p_value), ]
  
  ## Extract only the 'GO' terms from the results
  gprofiler_results_GOs <- gprofiler_results_reordered[grep('GO:', gprofiler_results_reordered$term_id), ]
  
  GOpval <- gprofiler_results_GOs[ , c("term_id", "p_value")]
  
  # Use REVIGO to visualize the results from the GO
  simMatrix <- calculateSimMatrix(GOpval$term_id,
                                  orgdb="org.Mm.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(GOpval$p_value), GOpval$term_id)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Mm.eg.db")
  file.name <- paste("treemap",i,sep = "_")
  file.name <- paste(file.name,"pdf",sep=".")
  pdf(file = paste(full.dir,file.name,sep = "/"), width = 8,height = 6)
  treemapPlot(reducedTerms)
  dev.off()
  
  mydf
}
