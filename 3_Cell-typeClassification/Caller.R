###Script for cell type calling
setwd("D:/users/fperez/NKI_TMAs_AF/R_tribus/")


for(sample.name in samples) {
  cat("Starting",sample.name,"...\n")
  t0 <- Sys.time()
  folder.name <- paste0("plots/", sample.name,"_gating")
  dir.create(folder.name, showWarnings = FALSE)
  
  # First call global cell types
  globalTypes <- cellTypeCaller(mydata[[sample.name]], global.gates, "GlobalCellType", folder.name = folder.name)
  firstGate <- merge(x=mydata[[sample.name]], y=globalTypes, all.x=T, by="CellId")
  
  
  # Then call subtypes of immune cells
  idx <- grep('Immune', firstGate$GlobalCellType)
  firstGate$GlobalCellType[idx] <- "Immune.cells"
  print(table(firstGate$GlobalCellType))
  
  sub.df <- firstGate[which(firstGate$GlobalCellType == "Immune.cells"),]
  if(nrow(sub.df)<100) {
    cat("cellType label has less than 100 immune cells.\n")
    next
  }else {
    cat(nrow(sub.df), "Immune cells passed to next gate.\n")
  }
  # Segond gate
  immuneTypes <- cellTypeCaller(sub.df, immune.gates, "ImmuneCellType", folder.name = folder.name)
  secondGate <- merge(x=firstGate, y=immuneTypes, all.x=T, by="CellId")
  print(table(immuneTypes$ImmuneCellType))
  
  sub.df <- secondGate[which(secondGate$ImmuneCellType == "CD8.T.cells"),]
  if(nrow(sub.df)<100) {
    cat("ImmuneCellType label has less than 100 CD8+ cells.\n")
    thirdGate <- secondGate
  }else {
    cat(nrow(sub.df), "CD8+ cells passed to next gate.\n")
    # Third gate
    cd8Types <- cellTypeCaller(sub.df, cd8.gates, "CD8TCellType", folder.name = folder.name)
    thirdGate <- merge(x=secondGate, y=cd8Types, all.x=T, by="CellId")
    print(table(cd8Types$CD8TCellType))
  }
  
  sub.df <- secondGate[which(secondGate$ImmuneCellType == "Macrophages"),]
  if(nrow(sub.df)<100) {
    cat("ImmuneCellType label has less than 100 Macrophages.\n")
    fourthGate <- thirdGate
  }else {
    cat(nrow(sub.df), "Macrophages passed to next gate.\n")
    # Third gate
    macsTypes <- cellTypeCaller(sub.df, macrophage.gates, "MacrophageType", folder.name = folder.name)
    fourthGate <- merge(x=thirdGate, y=macsTypes, all.x=T, by="CellId")
    print(table(macsTypes$MacrophageType))
  }
  
  sub.df <- secondGate[which(secondGate$ImmuneCellType == "CD4.T.cells"),]
  if(nrow(sub.df)<100) {
    cat("ImmuneCellType label has less than 100 CD4 T-cells\n")
    fifthGate <- fourthGate
  }else {
    cat(nrow(sub.df), "CD4 T-cells passed to next gate.\n")
    # Third gate
    cd4Types <- cellTypeCaller(sub.df, cd4.gates, "CD4Type", folder.name = folder.name)
    fifthGate <- merge(x=fourthGate, y=cd4Types, all.x=T, by="CellId")
    print(table(cd4Types$CD4Type))
  }
  
  print(table(fifthGate$ImmuneCellType)/sum(table(fifthGate$ImmuneCellType)))
  
  labeled[[sample.name]] <- fifthGate
  
  print.xy.plot(fifthGate, sample.name)
  # Save cellType column
  write.table(fifthGate, paste0("gated/",sample.name,"_gated.csv"), row.names = F)
  cat('Elapsed time',Sys.time() - t0,'seconds.\n')
  tryCatch({dev.off()},error=function(cond){return(NA)})
}
## Auxiliary function for sanity checks
print.xy.plot <- function(df, sample.name){
  df$cell.type <- df$GlobalCellType
  df$cell.type[which(df$GlobalCellType=="Immune.cells")] <- df$ImmuneCellType[which(df$GlobalCellType=="Immune.cells")]
  df$cell.type[ which(df$cell.type == "Other")] <- NA
  df$cell.type[ which(df$cell.type == "Unknown")] <- NA
  
  table(df$cell.type)
  ntypes <- length(unique(df$cell.type))
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(ntypes)
  
  p <- ggplot(df, aes(x=X_position, y=Y_position, color=cell.type)) +
    geom_point(size=1, stroke=0, shape='.') +
    scale_color_manual(values=mycolors, na.value = "grey50") +
    theme_bw() + ggtitle(sample.name) + coord_fixed(ratio = 1) +
    scale_y_reverse() +
    guides(colour = guide_legend(override.aes = list(size=8, shape=16)))
  gname <- paste0("XY_", sample.name, "scatter_by_immuneCellType.pdf")
  ggsave(plot=p, filename = gname, device="pdf")
}