### Additional functions for reproduction of analyses presented in Wilk, et al. Nature Medicine (2020) ###


add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


new_dotplot <- function(object = NULL, features = NULL, group.by = NULL, genes.on.x = TRUE, 
                        size.breaks.values = NULL, color.breaks.values = c(-3, -2, -1, 0, 1, 2, 3), shape.scale = 12, 
                        dend_x_var = "Average expression", dend_y_var = "Average expression",
                        cols.use = c("lightgrey", "blue"), scale.by = "radius", col.min = -2.5, col.max = 2.5,
                        dot.min = 0) {
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  
  data.features <- FetchData(object = object, vars = features)
  object[[group.by, drop = TRUE]]
  data.features$id <- object[[group.by, drop = TRUE]]
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if(genes.on.x){
    data.final <- data.frame(features.plot = data.plot$features.plot, id = data.plot$id)
    data.final <- cbind(data.final, data.plot[,c(2,5)])
  }
  else {
    data.final <- data.frame(id = data.plot$id, features.plot = data.plot$features.plot)
    data.final <- cbind(data.final, data.plot[,c(2,5)])
  }
  colnames(data.final)[3:4] <- c("Percent expressed", "Average expression")
  dot_plot(data.final, size_var = "Percent expressed", "Average expression", 
           dend_y_var = dend_y_var, dend_x_var = dend_x_var,
           dist = "euclidean", hclust_method = "ward.D2", x.lab.pos = "bottom",
           display_max_sizes = FALSE, size.breaks.values = size.breaks.values,
           shape.scale = shape.scale, color.breaks.values = color.breaks.values, cols.use = cols.use)
}


filter_quality <- function(.data) {
  final <- .data[-c(grep("^RPS", .data$gene), 
                    grep("^RPL", .data$gene), 
                    grep("^MT-", .data$gene),
                    grep("^MTR", .data$gene),
                    grep("MALAT1", .data$gene),
                    grep("^RNA18S5", .data$gene),
                    grep("^RNA28S5", .data$gene)),] 
  return(final)
}

plotPosCells <- function(seu, emat, features, grep = F) {
  if (grep) {
    pos <- emat[grep(paste(features,collapse="|"), rownames(emat)),]
    pos.cells <- colnames(pos[,colSums(pos) != 0])
  }
  else {
    if (length(features) >1) {
      pos <- emat[features,]
      pos.cells <- colnames(pos[,colSums(pos) != 0])
    }
    else {
      pos <- emat[features,] != 0
      pos.cells <- names(pos)[pos==T]
    }
    
  }
  
  DimPlot(seu, cells.highlight = pos.cells) + NoLegend()  
}


#non-numeric Abundance function
covid.seuAbundances <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"),
                                 meta.include = NULL, group_by = NULL, shape_by = NULL,
                                 custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                 pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                 ncol = 4, label = "p.signif", select.idents = NULL, 
                                 label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu@meta.data$cell.type, seu@meta.data[,"orig.ident"]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               "orig.ident"))  
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,"orig.ident"]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               "orig.ident"))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=uniques["seurat_clusters"]])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    
    
  }
  
} 


#numeric abundance function

covid.seuAbundances.num <- function (seu, by = c("orig.ident", "seurat_clusters", "cell.type"),
                                 meta.include = NULL, group_by = NULL, shape_by = NULL,
                                 custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                 pb = FALSE, correct = FALSE, comparisons = my_comparisons, 
                                 ncol = 4, label = "p.signif", select.idents = NULL, 
                                 label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)
  if (is.null(group_by)){
    group_by <- "null.group" 
  } 
  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 17, 15, 3, 7, 8)
    if ((length(unique((seu[[shape_by]]))[[1]])) > 6) {
      if (n > 18) {
        message(paste("At most 17 shapes are currently supported", 
                      "but", n, "are required. Setting 'shape_by' to NULL."))
        shape_by <- NULL
      }
      else {
        new <- setdiff(c(seq_len(16) - 1, 18), shapes)
        shapes <- c(shapes, new[seq_len(n - 6)])
      }
    }
  }
  if (by=="cell.type") {
    fq <- prop.table(table(seu@meta.data$cell.type, seu@meta.data[,"orig.ident"]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", 
                                                               "orig.ident"))  
  }
  else {
    fq <- prop.table(table(seu@meta.data$seurat_clusters, seu@meta.data[,"orig.ident"]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c("seurat_clusters", 
                                                               "orig.ident"))  
  }
  if (correct) {
    df <- df[grep("covid", df$orig.ident),]
    df$WBC <- mapvalues(df$orig.ident, from = covid_metadata.c$orig.ident, 
                        to = as.numeric(as.character(covid_metadata.c$WBC)))
    df$WBC <- as.numeric(as.character(df$WBC))
    df$freq <- df$freq*df$WBC
    df$WBC <- NULL
  }
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques<=uniques["seurat_clusters"]])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  if(pb){
    df <- df[df$cell.type %in% unique(covid_pb$cell.type),]
  }
  if(!is.null(select.idents)) {
    df <- df[df$cell.type %in% select.idents,]
  }
  if(correct) {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Cells (x1000/uL)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                              panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                  color = NA), strip.text = element_text(face = "bold"), 
                                                                                                              axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else {
    p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL, 
                                                                 y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                            panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                color = NA), strip.text = element_text(face = "bold"), 
                                                                                                            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), 
               position = "fill", stat = "identity") + scale_fill_manual("cell.type", 
                                                                         values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    if(is.null(custom_fill_colors)) {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()), 
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                            aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)),
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = 4) + guides(fill = FALSE) + geom_point(position = position_jitter(width = 0.25), 
                                                                                      aes_string(x = group_by, y = "freq", color = group_by.point, 
                                                                                                 shape = shape_by)) + scale_shape_manual(values = shapes) + theme(panel.grid.major = element_line(color = "grey", size = 0.25)) 
             + geom_line(aes_string(x = group_by, color = group_by.point, group = shape_by)))
    }
    else {
      switch(by, orig.ident = p + facet_wrap(group_by, scales = "free_x") + 
               geom_bar(aes_string(x = "orig.ident", fill = "factor(seurat_clusters)"), 
                        position = "fill", stat = "identity") + scale_fill_manual("seurat_clusters", 
                                                                                  values = cluster_cols) + scale_y_continuous(expand = c(0, 
                                                                                                                                         0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank()) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors), 
             
             
             
             seurat_clusters = p + facet_wrap("seurat_clusters", scales = "free_y", 
                                              ncol = 4) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by, 
                                                                                                         color = group_by, fill = group_by), position = position_dodge(), 
                                                                                              alpha = 0.25, outlier.color = NA) + geom_point(position = position_jitter(width = 0.25), 
                                                                                                                                             aes_string(x = group_by, y = "freq", color = group_by, 
                                                                                                                                                        shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors),
             
             
             
             
             cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                        ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                           alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                          aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                     shape = shape_by)) + scale_shape_manual(values = shapes) + 
               theme(panel.grid.major = element_line(color = "grey", 
                                                     size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label)
    }
    
    if(is.numeric(df[,group_by])){
      
      
      p + facet_wrap("cell.type", scales = "free_y", ncol = ncol) + 
        guides(fill = FALSE) +
        geom_point(size = pt.size, position = position_jitter(width = 0.25), 
                   aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)) +
        scale_shape_manual(values = shapes) + 
        theme(panel.grid.major = element_line(color = "grey", size = 0.25)) +
        scale_color_manual(values = custom_fill_colors) + 
        scale_fill_manual(values = custom_fill_colors) + 
        ggpubr::stat_cor(label.x = label.x)
    }  
    
  }
  
} 

writeAnnData <- function(emat = NULL, nmat = NULL, seu = NULL, idents = NULL, obs.vars = c("orig.ident", "seurat_clusters", "cell.type"), name = "seu", dir = NULL) {
  subset <- Seurat:::subset.Seurat(seu, idents = idents)
  seu.emat <- emat[,which(colnames(emat) %in% colnames(subset@assays$RNA)), drop=FALSE] 
  seu.nmat <- nmat[,which(colnames(nmat) %in% colnames(subset@assays$RNA)), drop=FALSE]
  seu.emat <- seu.emat[which(rownames(seu.emat) %in% rownames(seu.nmat)),, drop = FALSE]
  seu.nmat <- seu.nmat[which(rownames(seu.nmat) %in% rownames(seu.emat)),, drop = FALSE]
  seu.nmat <- seu.nmat[match(rownames(seu.emat), rownames(seu.nmat)),]
  if (!identical(colnames(seu.emat), colnames(seu.nmat), FALSE, FALSE, FALSE, FALSE)) {
    stop("Process failed: colnames not identical")
  }
  if(!identical(rownames(seu.emat), rownames(seu.nmat), FALSE, FALSE, FALSE, FALSE)) {
    stop("Process failed: rownames not identical")
  }
  
  #seu.obs <- as.data.frame(cbind(subset@meta.data$orig.ident, subset@meta.data$seurat_clusters, subset@meta.data$cell.type))
  seu.obs <- subset@meta.data[,colnames(subset@meta.data) %in% obs.vars]
  colnames(seu.obs)[grep("seurat_clusters", colnames(seu.obs))] <- "clusters"
  #colnames(seu.obs) <- obs.vars
  rownames(seu.obs) <- rownames(subset@meta.data)
  seu.obs <- seu.obs[match(colnames(seu.emat), rownames(seu.obs)),]
  seu.obsm <- as.data.frame(Embeddings(subset, reduction = "umap"))
  seu.obsm <- seu.obsm[match(colnames(seu.emat), rownames(seu.obsm)),]
  seu.obs_names <- as.data.frame(colnames(seu.emat))
  seu.var_names <- as.data.frame(rownames(seu.emat))
  
  writeMM(t(seu.emat), file = paste0(dir,name,".emat.mtx"))
  writeMM(t(seu.nmat), file = paste0(dir,name,".nmat.mtx"))
  write.csv(seu.obs, file =  paste0(dir,name,".obs.csv"), row.names = FALSE)
  write.csv(seu.obsm, file =  paste0(dir,name,".obsm.csv"), row.names = FALSE)
  write.table(seu.obs_names, file = paste0(dir,name,".obs_names.txt"), row.names = FALSE, col.names = FALSE)
  write.table(seu.var_names, file = paste0(dir,name,".var_names.txt"), row.names = FALSE, col.names = FALSE)
  message("Finished writing data")
}


ID.covid.markers <- function(seu = covid_combined.nc, idents = NULL, p.cutoff = 0.05){
  markers.list <- list()
  sub <- subset(covid_combined.nc, idents = idents)
  donors.2 <- unique(covid_combined.nc$Donor.full)[-grep("^H",
                                                         unique(covid_combined.nc$Donor.full))]
  for (k in 1:length(donors.2)) {
    try({
      markers <- subset(sub, cells = c(grep(donors.2[k],
                                            sub$Donor.full),
                                       grep("Healthy",
                                            sub$Status))) %>%
        FindMarkers(ident.1 = "COVID", group.by = "Status") %>% rownames_to_column(var = "gene") %>%
        add_column(Donor = donors.2[k])
      markers.list[[k]] <- markers[-c(grep("^RPS", markers$gene), 
                                      grep("^RPL", markers$gene), 
                                      grep("^MT-", markers$gene),
                                      grep("^MTR", markers$gene),
                                      grep("MALAT1", markers$gene),
                                      grep("^RNA18S5", markers$gene),
                                      grep("^RNA28S5", markers$gene)),]
    })
  }
  markers.list <- ldply(markers.list, data.frame)
  
  markers.list.mat <- markers.list[markers.list$p_val_adj<p.cutoff,]
  markers.list.mat <- markers.list.mat[,c("gene", "avg_logFC", "Donor")]
  
  markers.list.mat <- reshape2::dcast(markers.list.mat, formula = gene~Donor, value.var = "avg_logFC")
  markers.list.mat[is.na(markers.list.mat)] = 0
  markers.list.mat <- markers.list.mat %>% column_to_rownames(var = "gene")
  colnames(markers.list.mat) <- gsub("covid_", "", colnames(markers.list.mat))
  return(markers.list.mat)
}


covid.markers.heatmap <- function(markers.matrix = NULL, color = c("blue", "white", "red"), 
                                  paletteLength = 100, fontsize_row = 10, color.rows = T,
                                  title = "log(FC)", add.flag = T, de.cutoff = 4, top.n = NULL, 
                                  repel.degree = 0, legend = T, annotation_legend = T,
                                  cellwidth = NA, cellheight = NA,
                                  save = F, width = 7, height = 11, file = "~/Downloads/p.pdf") {
  
  paletteLength = 100
  if(sum(markers.matrix>=0)==dim(markers.matrix)[1]*dim(markers.matrix)[2]) {
    color = color[2:3]
  }
  myColor = colorRampPalette(color)(paletteLength)
  myBreaks <- unique(c(seq(min(markers.matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(max(markers.matrix)/paletteLength, max(markers.matrix),
                           length.out=floor(paletteLength/2))))
  annotation_colors = list(
    Ventilated = c(ARDS="red3", NonVent=RColorBrewer::brewer.pal(9, "Oranges")[4]))
  p <- pheatmap(markers.matrix, color = myColor, breaks = myBreaks, 
                heatmap_legend_param = list(title = title), angle_col = "90", 
                annotation_col = row_annotation, annotation_colors = annotation_colors, legend = legend,
                annotation_legend = annotation_legend, fontsize_row = fontsize_row,
                cellwidth = cellwidth, cellheight = cellheight)
  if(color.rows){
    markers.matrix <- markers.matrix[match(p$gtable$grobs[[5]]$label,rownames(markers.matrix)),]
    p$gtable$grobs[[5]]$gp=gpar(col=ifelse((rowSums(markers.matrix))>0, "red", "blue"), fontsize = fontsize_row)
  }
  
  if(add.flag){
    if(save){
      pdf(file = file, width = width, height = height)
    }
    if(!is.null(top.n)) {
      kept.labels <- names((abs(rowSums(markers.matrix)) %>% sort(decreasing = T))[1:top.n])
    }
    else {
      kept.labels = names(rowSums(markers.matrix !=0)[rowSums(markers.matrix !=0)>=de.cutoff])
    }
    add.flag(p, kept.labels = kept.labels,
             repel.degree = repel.degree)
    if(save){
      dev.off()
    }
  }
  else{
    p
  }
}



build.cp.matrix <- function(compartment = NULL, by = "z", p.cutoff = 0.05, z.cutoff = T, top.n = 100, path = "/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/Aaron/Data/scRNA-seq/COVID/analyses/IPA/results/") {
  files <- list.files(path)
  files_import <- grep(paste0(compartment,".txt"),files,value = T)
  donor_names <- gsub(paste0("_",compartment,".txt"),"",files_import)
  donor_names <- gsub("covid_","",donor_names)
  datalist = lapply(files_import, function(x)read_delim(paste0("/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/Aaron/Data/scRNA-seq/COVID/analyses/IPA/results/",x), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1))
  datalist <- mapply(cbind, datalist, "Donor"=donor_names,SIMPLIFY = F)
  datalist.df <- ldply(datalist,data.frame)
  datalist.df$X6 <- NULL
  colnames(datalist.df) <- mapvalues(colnames(datalist.df), from = c("Ingenuity.Canonical.Pathways",
                                                                     "X.log.p.value.",
                                                                     "Ratio",
                                                                     "z.score",
                                                                     "Molecules",
                                                                     "Donor"),
                                     to = c("Pathway",
                                            "p",
                                            "Ratio",
                                            "z",
                                            "Molecules",
                                            "Donor"))
  datalist.final <- datalist.df
  if(!is.null(p.cutoff)){
    p.converted = -log(0.05,10)
    datalist.final <- datalist.final[datalist.final$p>p.converted,]
  }
  if(z.cutoff){
    datalist.final <- datalist.final[!is.na(datalist.final$z),]
    datalist.final <- datalist.final[!datalist.final$z==0,]
  }
  datalist.final[is.na(datalist.final$z),"z"] <- 0
  datalist.mat <- reshape2::dcast(data = datalist.final, formula = Pathway~Donor,fun.aggregate = sum, value.var = by) %>% column_to_rownames(var = "Pathway")
  if(by=="p" && z.cutoff){
    datalist.z <- as.matrix(reshape2::dcast(data = datalist.final, formula = Pathway~Donor,fun.aggregate = sum, value.var = "z") %>% column_to_rownames(var = "Pathway"))
    datalist.mat <- as.matrix(datalist.mat)
    datalist.mat[datalist.z<0] <- -datalist.mat
    datalist.mat <- as.data.frame(datalist.mat)
  }
  
  if(!is.null(top.n)) {
    top <- names((abs(rowSums(datalist.mat)) %>% sort(decreasing = T))[1:top.n])
    datalist.mat <- datalist.mat[rownames(datalist.mat) %in% top,]
  }
  colnames(datalist.mat) <- mapvalues(colnames(datalist.mat), from = c("555_1",
                                                                       "555_2",
                                                                       "556",
                                                                       "557",
                                                                       "558",
                                                                       "559",
                                                                       "560",
                                                                       "561"),
                                      to = c("C1 A",
                                             "C1 B",
                                             "C2",
                                             "C3",
                                             "C4",
                                             "C5",
                                             "C6",
                                             "C7"))
  return(datalist.mat)
}



build.ur.matrix <- function(compartment = NULL, by = "z", p.cutoff = 0.05, z.cutoff = T, top.n = 100, path = "/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/Aaron/Data/scRNA-seq/COVID/analyses/IPA/results/") {
  files <- list.files(path)
  files_import <- grep(paste0(compartment,".ur.txt"),files,value = T)
  donor_names <- gsub(paste0("_",compartment,".ur.txt"),"",files_import)
  donor_names <- gsub("covid_","",donor_names)
  datalist = lapply(files_import, function(x)read_delim(paste0("/Volumes/GoogleDrive/My Drive/Blish Lab/00 - All Server Data and Folders/Aaron/Data/scRNA-seq/COVID/analyses/IPA/results/",x), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1))
  datalist <- mapply(cbind, datalist, "Donor"=donor_names,SIMPLIFY = F)
  datalist.df <- ldply(datalist,data.frame)
  colnames(datalist.df) <- mapvalues(colnames(datalist.df), from = c("Upstream.Regulator",
                                                                     "Expr.Log.Ratio",
                                                                     "Molecule.Type",
                                                                     "Predicted.Activation.State",
                                                                     "Activation.z.score",
                                                                     "p.value.of.overlap",
                                                                     "Target.Molecules.in.Dataset",
                                                                     "Mechanistic.Network",
                                                                     "Donor",
                                                                     "Flags"),
                                     to = c("Regulator",
                                            "Expr",
                                            "Type",
                                            "State",
                                            "z",
                                            "p",
                                            "Molecules",
                                            "Network",
                                            "Donor",
                                            "Flags"))
  datalist.final <- datalist.df
  if(!is.null(p.cutoff)){
    datalist.final <- datalist.final[datalist.final$p<p.cutoff,]
  }
  if(z.cutoff){
    datalist.final <- datalist.final[!is.na(datalist.final$z),]
    datalist.final <- datalist.final[!datalist.final$z==0,]
  }
  datalist.final[is.na(datalist.final$z),"z"] <- 0
  datalist.mat <- reshape2::dcast(data = datalist.final, formula = Regulator~Donor,fun.aggregate = sum, value.var = by) %>% column_to_rownames(var = "Regulator")
  if(by=="p") {
    datalist.mat <- -log(datalist.mat,10)
  }
  if(!is.null(top.n)) {
    top <- names((abs(rowSums(datalist.mat)) %>% sort(decreasing = T))[1:top.n])
    datalist.mat <- datalist.mat[rownames(datalist.mat) %in% top,]
  }
  colnames(datalist.mat) <- mapvalues(colnames(datalist.mat), from = c("555_1",
                                                                       "555_2",
                                                                       "556",
                                                                       "557",
                                                                       "558",
                                                                       "559",
                                                                       "560",
                                                                       "561"),
                                      to = c("C1 A",
                                             "C1 B",
                                             "C2",
                                             "C3",
                                             "C4",
                                             "C5",
                                             "C6",
                                             "C7"))
  return(datalist.mat)
}

new.markers.heatmap <- function(markers.matrix = NULL, fontsize_row = 10, color.rows = T,
                                title = "log(FC)", add.flag = T, de.cutoff = 4, 
                                top.n = NULL, 
                                repel.degree = 0, legend = T, annotation_legend = T,
                                cellwidth = NA, cellheight = NA,
                                save = T, width = 7, height = 11, 
                                file = "~/Downloads/p.pdf") {
  
  ta = columnAnnotation(Age = row_annotation$Age, 
                        DPS = row_annotation$DPS,
                        DTF = row_annotation$DTF,
                        Admission = row_annotation$Admission,
                        ARDS = row_annotation$Ventilated,
                        col = list(ARDS = c("No" = "orange", "Yes" = "red"),
                                   Admission = c("Floor" = "blue", "ICU" = "magenta")))
  if(add.flag){
    if(!is.null(top.n)) {
      top.labels <- names((abs(rowSums(markers.matrix)) %>% sort(decreasing = T))[1:top.n])
      kept.labels <- rownames(markers.matrix)[rownames(markers.matrix) %in% top.labels]
    }
    else {
      kept.labels = names(rowSums(markers.matrix !=0)[rowSums(markers.matrix !=0)>=de.cutoff])
    }
    t = 1:nrow(markers.matrix)
    
    ha = rowAnnotation(foo = anno_mark(at = t[rownames(markers.matrix) %in% kept.labels], 
                                       labels = kept.labels,
                                       labels_gp = gpar(col =
                                                          ifelse((rowSums(markers.matrix[rownames(markers.matrix) %in% kept.labels,]))>0,
                                                                 "red", "blue"))))
  }
  
  
  p <- Heatmap(markers.matrix, right_annotation = ha, 
               top_annotation = ta, row_names_gp = gpar(fontsize = 0))
  if(save){
    pdf(file = file, width = width, height = height)
    p
    #dev.off()
  }
  else {
    p
  }
}

processNewSeurat <- function(parent.object, idents = NULL, cells = NULL) {
  if (!is.null(idents)) {
    seu <- subset(parent.object, idents = idents)
  }
  if (!is.null(cells)) {
    seu <- subset(parent.object, cells = cells)
  }
  message("Running PCA")
  seu <- RunPCA(seu, verbose = FALSE)
  message("Running UMAP")
  seu <- RunUMAP(seu, dims = 1:50, verbose = FALSE)
  message("Clustering")
  seu <- FindNeighbors(seu, dims = 1:50, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
  return(seu)
}



myFeaturePlot <- function(object = NULL, features = features, ncol = NULL, save = T, save.as = "png", height = 7, width = 7, cols = c("lightgrey", "blue")) {
  plots <- lapply(features, function(x) {FeaturePlot(object,features = x, cols = cols) +  labs(x = "UMAP1", y = "UMAP2") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), axis.line = element_blank())})
  CombinePlots(plots = plots, ncol = ncol)
  if(save) {
    if(save.as=="pdf"){
      ggsave("p.pdf", path = "~/Downloads/", height = height, width = width)
    }
    if(save.as=="png"){
      ggsave("p.png", path = "~/Downloads/", height = height, width = width)
    }
  }
}

myHighlightCells <- function(object = NULL, idents = NULL, group_by = "seurat_clusters", ncol = NULL, save = T, height = 7, width = 7, col = "black") {
  cells.highlight=list()
  for(i in 1:length(idents)){
    cells = colnames(object)[grep(idents[i], object@meta.data[,group_by])]
    cells.highlight[[i]]=cells
  }
  
  plots <- lapply(cells.highlight, function(x) {
    DimPlot(object,cells.highlight = x, cols.highlight = col) + 
      NoLegend() + labs(x = "UMAP1", y = "UMAP2") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  })
  CombinePlots(plots = plots, ncol = ncol)
  if(save) {
    ggsave("p.pdf", path = "~/Downloads/", height = height, width = width)
  }
}


constructConsensus <- function(markers.matrix, donors.use = 1:8, save = T, height = 5, width = 7) {
  markers.matrix.m <- markers.matrix %>% rownames_to_column(var = "gene")
  markers.matrix.m <- reshape2::melt(markers.matrix.m)
  order.markers <- rownames(markers.matrix[order(rowSums(-markers.matrix)),])
  markers.sum <- as.data.frame(ifelse(markers.matrix != 0, 1, 0))
  markers.sum$total <- rowSums(markers.sum)
  keep <- rownames(markers.sum[markers.sum$total>3,])
  markers.matrix.m <- markers.matrix.m[markers.matrix.m$gene %in% keep,]
  ggplot(markers.matrix.m, aes(x = factor(gene, level = order.markers), 
                               y = value, fill = variable)) +
    geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() + 
    labs(x = "", y = "Cumulative average log(fold-change)", fill = "Sample") +
    scale_fill_manual(values = custom_fill_colors[donors.use]) + 
    ggpubr::rotate_x_text()
  if(save){
    ggsave("p.pdf", path = "~/Downloads/", height = height, width = width)
  }
}







