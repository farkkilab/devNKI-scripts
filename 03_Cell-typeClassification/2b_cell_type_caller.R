cellTypeCaller <- function(df, gates, gates.class="", folder.name=folder.name, grid.size.xdim=15, hierarchical.trees=TRUE, scaling=TRUE) {
  marker_cols <- unique(unlist(gates))
  mat <- as.matrix(df[, marker_cols])
  if (scaling){
    print ("Starting to scale data")
    mat <- scale(mat)
    mat <- apply(mat, 2, min_max_norm)
  }
  data_FlowSOM <- flowCore::flowFrame(mat)
  
  n.cell.types <- length(gates)
  # set seed for reproducibility
  set.seed(42)
  
  ## run FlowSOM ###
  out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
  out <- FlowSOM::BuildSOM(out, colsToUse = marker_cols, silent=T, xdim=grid.size.xdim, ydim=grid.size.xdim)
  out <- FlowSOM::BuildMST(out, silent=T)
  # save hierarchical trees colored by markers
  if (hierarchical.trees){
    pdf(paste0(folder.name, "/",gates.class,"_marker_expr_by_node.pdf"))
    for(chan in marker_cols){
      PlotMarker(out, chan)
    }
    tryCatch({dev.off()},error=function(cond){return(NA)})
  }
  
  ### Gate scores for each node ##
  gate.scores.by.node <- list()
  tryCatch({dev.off()},error=function(cond){return(NA)})
  
  pdf(paste0(folder.name, "/",gates.class,"_gate_score_by_node.pdf"))
  for(gate.name in names(gates)){
    # Transform user list into valid query format
    query <- c(rep("high", length(gates[[gate.name]]$Pos)), rep("low", length(gates[[gate.name]]$Neg)))
    names(query) <- c(gates[[gate.name]]$Pos, gates[[gate.name]]$Neg)
    # Run query for each node
    query_res = QueryStarPlot(out, query, plot = FALSE, main=gate.name)
    res <-query_res$nodeScores ## Changed in new
    res[!query_res$selected] <- min(res)
    gate.scores.by.node[[gate.name]] <- res
  }
  tryCatch({dev.off()},error=function(cond){return(NA)})
  #Concatenate celltype scores in dataframe
  gate.scores <- do.call(cbind, gate.scores.by.node)
  nodeID.per.cell <- out$map$mapping[,1]
  annotation.scores <- data.frame(gate.scores, row.names = 1:nrow(gate.scores))

  # Median expression of the gating channels to annotate the heatmap
  medians.by.node <- aggregate(mat, by = list(out$map$mapping[,1]), median)
  
  ann_row <- data.frame(scale(medians.by.node[,-1]), row.names = 1:nrow(medians.by.node))
  if (ncol(ann_row) > 15 ) {
    ann_row <- ann_row[, order(colnames(ann_row))[1:15]]
  }
  tryCatch({dev.off()},error=function(cond){return(NA)})
  if (ncol(annotation.scores) == 1) annotation.scores$Unknown <- 1-annotation.scores[,1]
  p <- pheatmap(annotation.scores,
                show_rownames=T,show_colnames=T,cluster_rows=T,cluster_cols=T,method="ward.D2",
                cellheight=4,cellwidth=10,fontsize_row=4,fontsize_col=12,
                color=colorRampPalette(c("blue","white","red"),interpolate="linear")(300),
                border_color="white", scale='none',
                annotation_row = ann_row, cutree_rows =max(nodeID.per.cell), annotation_legend = FALSE,
                filename =paste0(folder.name, "/",gates.class,"_normalized_scores_heatmap.pdf"))

  # Set labels based on clustering of scores and highest scored gate on each node
  class.labels <- apply(annotation.scores, 1, function(x) {
    if(max(x) <= 0.5)
      "Unknown"
    else
      colnames(annotation.scores)[which.max(x)]
  })
  
  tryCatch({dev.off()},error=function(cond){return(NA)})
  # Asign cell type to each cell based on the node and the gate for that node
  nodeID.per.cell <- out$map$mapping[,1]
  df[,gates.class] <- class.labels[nodeID.per.cell]
  
  # Plot resulting tree and compare manually to the marker trees
  pdf(paste0(folder.name,  "/",gates.class,"_cellType_tree.pdf"))
  PlotStars(out, backgroundValues=class.labels)
  tryCatch({dev.off()},error=function(cond){return(NA)})
  
  res1 <- df[,c('CellId',gates.class)]
  rest.return <- list(res1, annotation.scores, nodeID.per.cell)
  return(rest.return)
}
