######################Function to make complex heatmaps##################
heatmaps <- function(data.cell, method=NULL){
  # cells_percentages <- data.cell %>% group_by(imageid, rcn_id) %>% summarise(n = n()) %>% mutate(proportion = n * 100 / sum(n))
  # cells_percentages <- as.data.frame(cells_percentages)
  # cells_percentages <- cells_percentages[,-3]
  # 
  # cells_percentages.pt <- merge(cells_percentages, clindata.shrd, by="imageid")[,c(1:3,6)]
  
  # cells_percentages.pt <- cells_percentages.pt %>%
  #   group_by(patient, rcn_id) %>%
  #   summarise(m.proportion = mean(proportion))
  # cells_percentages <- as.data.frame(cells_percentages.pt)
  
  cells_percentages.pt <- data.cell %>%
    group_by(patient, rcn_id) %>%
    summarise(n = n()) %>%
    mutate(m.proportion = n * 100 / sum(n))
  cells_percentages <- as.data.frame(cells_percentages.pt)
  cells_percentages <- cells_percentages[,-3]
  
  cells_percentages.m <- cells_percentages %>% pivot_wider(names_from = rcn_id, values_from = m.proportion)
  
  cells_percentages.m <- as.data.frame(cells_percentages.m)
  row.names(cells_percentages.m) <- cells_percentages.m[,1]
  cells_percentages.m <- cells_percentages.m[,-1]
  cells_percentages.m[is.na(cells_percentages.m)] <- 0
  
  #Trimming outliers to the 99 percentile
  cells_percentages.m <- apply(cells_percentages.m,2,function(x){
    max.to.trim <- quantile(x, 0.99)
    x[which(x > max.to.trim)] <- max.to.trim
    return(x)
  })
  cells_percentages.scaled <- scale(cells_percentages.m)
  
  #Getting molecular profile of patients using same order as matrix
  Profiles.m <- sapply(row.names(cells_percentages.m), function(x) {
    clindata.shrd[which(clindata.shrd$patient == x)[1],"Molecular.profile2"] })
  Profiles.m <- factor(Profiles.m, levels=c("HRD", "BRCAmut/met","BRCAness-sHRD+","BRCAness+sHRD-","CCNE1amp","HRP","Other"))
  
  #Getting therapy sequence of patients using same order as matrix
  therapy <- sapply(row.names(cells_percentages.m), function(x) {
    clin.tcycif[which(clin.tcycif$patient == x)[1],"therapy_sequence"] })
  therapy[which(therapy == "PDS followed by NACT")] <- "PDS"
  therapy[which(therapy == "Primairy debulking")] <- "PDS"
  therapy[which(therapy == "Only debulking")] <- "PDS"
  therapy[which(therapy == "NACT followed by re-debulking")] <- "NACT"
  
  #colours <- list(Therapy=c("NACT"="black", "PDS"="grey90"),
  #                M.profile=c("BRCAmut/met" = "#ff0000", "HRD" = "#960018", "CCNE1amp"="darkblue", "HRP"="#1E90FF", "BRCAness+sHRD-"="pink", "BRCAness-sHRD+"="#BA3405",  "Other"="grey60"))
  
  colours <- list(Therapy=c("NACT"="black", "PDS"="grey90"),
                  M.profile=c("BRCAmut/met" = "#ff0000", "HRD" = "#960018", "CCNE1amp"="darkblue", "HRP"="#1E90FF", "Other"="grey60"))
  
  ha = rowAnnotation(Therapy=therapy, M.profile=Profiles.m, col =colours)
  
  set.seed(123)
  col_c <- colorRamp2(c(-1.5, 0, 2, 4), c("#0000FF","grey95", "#F9AAAA", "#EE0000"))
  hmap1 <- Heatmap(cells_percentages.scaled, name="Subtype %\nper patient (Z-score)",
                   row_names_gp = gpar(fontsize = 2), column_names_gp = gpar(fontsize = 7), row_title = NULL, column_title = NULL,
                   row_names_side = "left", row_dend_side = "right",
                   col=col_c, #row_km = 3, column_km = 4,
                   right_annotation = ha, column_names_rot = 60)
  png(paste0(out.folder.name, "RCN_all-heatmap_", method, ".png"), res=300, width=5, height=7, units = "in")
  draw(hmap1)
  dev.off()
  
  therapies.vals <- c("NACT","PDS")
  for (t in therapies.vals){
    cells_percentages.sel.scaled <- cells_percentages.scaled[which(therapy == t),]
    
    ha = rowAnnotation(M.profile=Profiles.m[which(therapy == t)], col=colours)
    
    set.seed(321)
    col_c <- colorRamp2(c(-1.5, 0, 2, 4), c("#0000FF","grey95", "#F9AAAA", "#EE0000"))
    hmap1 <- Heatmap(cells_percentages.sel.scaled, name="Subtype %\nper patient (Z-score)",
                     row_names_gp = gpar(fontsize = 3), column_names_gp = gpar(fontsize = 7), column_title = t, row_title = NULL,
                     row_names_side = "left", row_dend_side = "right",
                     col=col_c, #row_km = 3, column_km = 3,
                     right_annotation = ha, column_names_rot = 60)
    
    png(paste0(out.folder.name, t, "_RCN_heatmap_", method, ".png"), res=300, width=5, height=7, units = "in")
    draw(hmap1)
    dev.off()
  }
}


###################Function for calculating diversity and abundance for each RCN##################
heatmap.diversity <- function(input.cells, method="NA"){
  cells_percentages <- input.cells %>% 
    group_by(rcn_id, GlobalCellType) %>% 
    summarise(n = n()) %>% mutate(count = n)
  cells_percentages <- as.data.frame(cells_percentages)
  cells_percentages <- cells_percentages[,-3]
  
  cells_percentages.m <- cells_percentages %>% 
    pivot_wider(names_from = GlobalCellType, values_from = count)
  
  cells_percentages.m <- as.data.frame(cells_percentages.m)
  row.names(cells_percentages.m) <- cells_percentages.m[,1]
  cells_percentages.m <- cells_percentages.m[,-1]
  cells_percentages.m[is.na(cells_percentages.m)] <- 0
  
  H <- rev(diversity(cells_percentages.m))
  cluster_abundances <- as.vector(rev(table(input.cells$rcn_id)))
  
  l <- length(H)
  
  col_fun = colorRamp2(c(min(H), max(H)), c("pink", "darkred"))
  
  Rowann <- rowAnnotation(Diversity=H, col = list(Diversity = col_fun),
                          cells = anno_barplot(cluster_abundances))
  
  ht <- Heatmap(matrix(rnorm(l * 3), l), left_annotation = Rowann, cluster_rows = FALSE)
  
  svg(paste0(out.folder.name,"Diversity_RCN_heatmap_", method, ".svg"), width = 6, height = 5.4)
  draw(ht)
  dev.off()
}

##############Just calculate the proportions of cells by global.cell.type##################
counts.by.GlobalCellType <- function(input.cells){
  cells_percentages <- input.cells %>% 
    group_by(rcn_id, GlobalCellType) %>% 
    summarise(n = n()) %>% mutate(count = n / sum(n))
    #summarise(n = n()) %>% mutate(count = n)
  cells_percentages <- as.data.frame(cells_percentages)
  cells_percentages <- cells_percentages[,-3]
  
  cells_percentages.m <- cells_percentages %>% 
    pivot_wider(names_from = GlobalCellType, values_from = count)
  
  cells_percentages.m <- as.data.frame(cells_percentages.m)
  row.names(cells_percentages.m) <- cells_percentages.m[,1]
  cells_percentages.m <- cells_percentages.m[,-1]
  cells_percentages.m[is.na(cells_percentages.m)] <- 0
  return(cells_percentages.m)
}


##############Barplot for the RCN composition##################
barplot.rcn <- function(input.cells, method="NA"){
  cells_percentages <- input.cells %>%
    group_by(rcn_id, GlobalCellType) %>%
    summarise(n = n()) %>% mutate(proportion = n * 100 / sum(n))
  cells_percentages <- as.data.frame(cells_percentages)
  cells_percentages <- cells_percentages[,-3]
  
  #GC_order <- c("CancerC1","CancerC2","CancerC3","CancerC4","CancerC5","CancerC6","CancerC7","CancerC8","B.cells","T.regs","CD4.T.cells","CD8.T.cells","CD68.MP","CD15.MY","CD11c.MY","CD163.MP","Other.MY","Other.immune","StromalC1","StromalC2","StromalC3","StromalC4","StromalC5","StromalC6","StromalC7","StromalC8")
  GC_order <- c("Cancer","B.cells","T.regs","CD4.T.cells","CD8.T.cells",
                "CD68.MP","CD15.MY","CD11c.MY","CD163.MP","CD207.MY","Other.MY",
                "Other.immune", "Stromal","Others")
  
  cells_percentages$GlobalCellType <- factor(cells_percentages$GlobalCellType, levels=GC_order)
  
  p <- ggplot(cells_percentages, aes(x=rcn_id,y=proportion, fill=GlobalCellType)) +
    geom_bar(stat="identity", color="grey20") + theme_classic() +
    theme(axis.text.x=element_text(size=rel(1.2))) + coord_flip() + 
    xlab("") + ylab("Proportion") +
    scale_fill_manual(values = as.character(paletteer_c("grDevices::rainbow", length(GC_order))))
  print(paste0("Saving file in:", paste0(out.folder.name,"RCN_barplot_", method,".png")))
  ggsave(p, file=paste0(out.folder.name,"RCN_barplot_", method,".png"), width = 14, height = 12, units = "cm")
  ggsave(p, file=paste0(out.folder.name,"RCN_barplot_", method,".svg"), width = 14, height = 12, units = "cm")
}
