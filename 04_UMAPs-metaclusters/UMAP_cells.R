
############### DATA ############################################################################################

#For Figure 1
library(viridis)
library(ggplot2)
library(uwot)
library(RColorBrewer)
library(gridExtra)
library(grid) 

source("D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/4_UMAPs/Z_trimming.R")

#Number of random cells for UMAP
N <- 50000
#ID of cells for iterations
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(8)
#subtypes.labels <- c("Tumor", "Stromal.cells", "Lymphoid1", "Lymphoid2", "Myeloid.cells",  "LostCell", "LowDapi", "Other")
subtypes.labels <- c("Tumor", "Stromal.cells", "Lymphoid1", "Lymphoid2", "Myeloid.cells")

source("D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/3_Cell-typeClassification/gates_NKI.R")

samples <- c("TMA_18_810","TMA_33_576","TMA_45_312")
names.dataset <- c("0_expansion", "1px_expansion", "2px_expansion", "3px_expansion")

for (j in 1:length(samples)){
  print(paste0("Starting to analyze sample:", samples[j]))
  
  print("Reading input files")
  celltypes0 <- read.csv(paste0("D:/users/fperez/NKI_TMAs_AF/",samples[j], "/Tribus_Nuclei-segment/", samples[j], "_cellTypes.csv"))
  celltypes1 <- read.csv(paste0("D:/users/fperez/NKI_TMAs_AF/" ,samples[j], "/Tribus_Cell-segment1/", samples[j], "_cellTypes.csv"))
  celltypes2 <- read.csv(paste0("D:/users/fperez/NKI_TMAs_AF/" ,samples[j],"/Tribus_Cell-segment2/", samples[j], "_cellTypes.csv"))
  celltypes3 <- read.csv(paste0("D:/users/fperez/NKI_TMAs_AF/" ,samples[j], "/Tribus_Cell-segment3/", samples[j], "_cellTypes.csv"))
  
  print("Starting preprocesing")
  celltypes0 <- celltypes0[!(celltypes0$GlobalCellType %in% c("LostCell", "LowDapi", "Other")),]
  celltypes1 <- celltypes1[!(celltypes1$GlobalCellType %in% c("LostCell", "LowDapi", "Other")),]
  celltypes2 <- celltypes2[!(celltypes2$GlobalCellType %in% c("LostCell", "LowDapi", "Other")),]
  celltypes3 <- celltypes3[!(celltypes3$GlobalCellType %in% c("LostCell", "LowDapi", "Other")),]
  
  cell.annotations <- list(celltypes0, celltypes1,  celltypes2, celltypes3)
  
  gate.cols <- colnames(celltypes0)[which(colnames(celltypes0) %in% unique(unlist(global.gates)))]
  random.cells <- sort(sample(1:nrow(celltypes0), N))
  
  ################   GLOBAL CELL TYPE UMAPS ################################
  
  print("Performing UMAPs")
  pltList <- list()
  for (i in c(1:4)){
    #Selecting subset of channels and cells
    all_celltypes <- cell.annotations[[i]]
    df <- all_celltypes[, c("CoreId",gate.cols ,"GlobalCellType")]
    colnames(df)[15] <- "Subtype"
    df.subset  <- df[random.cells,]
    df.subset[,c(2:14)] <- log2(df.subset[,c(2:14)])
    #data.truncated <- z.trimming(df.subset[,c(2:14)])
    #df.subset[,c(2:14)] <- data.truncated[[1]]
    
    df.subset$Subtype <- factor(df.subset$Subtype, levels=subtypes.labels)
    cell_type <- df.subset$Subtype
    
    umap_s = uwot::umap(df.subset[,c(2:14)],n_neighbors = 80, scale=T ,spread=1, min_dist = 0.5,n_epochs = 60)
    uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
    
    if (i == 1){
      pltChannels <- list()
      pltChannels <- lapply(2:14, function(x){
        markers <- df.subset[, x]
        p = ggplot(data=uplot,aes(x,y, color=markers)) + geom_point(size=0.3, stroke=0, alpha=0.7)+ 
          guides(color = guide_legend(override.aes = list(size=5))) + 
          theme_bw() + coord_fixed(ratio=1)+ 
          xlab("umap1") + ylab("umap2") + theme(legend.position="none") + theme(aspect.ratio=1)
        p <- p + ggtitle(colnames(df.subset)[x])
        print(p)
      })
      ggsave(file=paste0("D:/users/fperez/NKI_TMAs_AF/",samples[j],"_UMAPs_channels_untruncated.pdf"), arrangeGrob(grobs = pltChannels, ncol = 4), height = 14, width = 16)
    }
    
    p = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
      geom_density_2d(color="black") +
      scale_fill_viridis(discrete = TRUE) + 
      theme(legend.title = element_blank()) + guides(color = guide_legend(override.aes = list(size=5))) + 
      scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
      guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="bottom", byrow=F)) + 
      xlab("umap1") + ylab("umap2") + #ylim(-10, 10) + 
      theme(legend.position="none") + theme(aspect.ratio=1)
    p <- p + ggtitle(names.dataset[i])
    print(p)
    pltList[[i]] <- print(p)
  }
  ggsave(file=paste0("D:/users/fperez/NKI_TMAs_AF/",samples[j],"_UMAPs_untruncated.pdf"), arrangeGrob(grobs = pltList, ncol = 4), height = 7, width = 17)
}


####################### COLORED BY PATIENT OF ORIGIN #########################################################
Patient <- as.factor(df$Sample)
n <- 44
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycolors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p1 = ggplot(uplot,aes(x,y, color=Patient))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  theme_bw() + coord_fixed(ratio=1)+scale_color_manual(values=mycolors) +
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=3, label.position="bottom", byrow=F)) + 
  xlab("umap1") + ylab("umap2") + theme(legend.position="none", aspect.ratio=1) + xlim(-10, 10) + ylim(-10, 10)
print(p1)


##############################################################################################################
##################### IMMUNE CELL UMAPS ######################################################################

#all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")

all_celltypes_24092020 <- all_celltypes_immune
all_celltypes_24092020$Subtype <- as.character(all_celltypes_24092020$Subtype)
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD11c+CD163+IBA1+Macrophages"), "Subtype"] <- "IBA1+CD163+CD11c+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD11c+IBA1+Macrophages"), "Subtype"] <- "IBA1+CD11c+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD163+IBA1+Macrophages"), "Subtype"] <- "IBA1+CD163+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD4"), "Subtype"] <- "CD4+ Effector T-cells"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD8"), "Subtype"] <- "CD8+ T-cells"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD163+Macrophages"), "Subtype"] <- "CD163+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "IBA1+Macrophages"), "Subtype"] <- "IBA1+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "CD11c+APC"), "Subtype"] <- "CD11c+APC"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "FOXP3+CD4+Tregs"), "Subtype"] <- "FOXP3+CD4+ T-regs"


#all_cell_types_16092020_immune <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Immune"),]

df <- all_cell_types_16092020_immune[, c("Sample","CD163","CD20","CD4","CD3d","CD8a",
                                         "FOXP3","IBA1","CD11c", "CD1c", "Subtype")]
df$Subtype <- as.factor(df$Subtype)


cell_type <- df$Subtype

umap_s = uwot::umap(df[,c(2:10)],n_neighbors = 80, scale=T ,spread=2.5,n_epochs = 60, min_dist = 0.2)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
mycolors <- colorRampPalette(brewer.pal(10, "Paired"))(10)
p2 = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.7, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="right", byrow=F, reverse=TRUE)) + 
  xlab("umap1") + ylab("umap2")+ theme(legend.position="none", aspect.ratio=1)
print(p2)

############################ COLORED BY MARKER ################################################################

for (i in c(2:10)){
  Expression <- df[, i]
  p = ggplot(uplot,aes(x,y, color=Expression))+ geom_point(size=0.5, stroke=0, alpha=0.7)+  
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
}


################################################################################################################
####################### STROMAL CELL UMAPS #####################################################################

#all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")

all_celltypes_24092020 <- all_celltypes_stroma
all_celltypes_24092020$Subtype <- as.character(all_celltypes_24092020$Subtype)
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "High-PDL1"), "Subtype"] <- "Functional stroma"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Non-proliferative_Stroma"), "Subtype"] <- "Non-proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Low_eccentricity_medium_vimentin"), "Subtype"] <- "Low eccentricity"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Proliferative_Stroma"), "Subtype"] <- "Proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "High-proliferative_Stroma"), "Subtype"] <- "High-proliferative Stroma"

#all_cell_types_13102020_stromal <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Stroma_Endothelia"),]

df <- all_celltypes_24092020[, c("Sample","cCasp3", "pSTAT1", "Ki67",
                                 "PDL1","vimentin","CD31","P21","Eccentricity", "Subtype")]
df$Subtype <- as.factor(df$Subtype)


cell_type <- df$Subtype
umap_s = uwot::umap(df[,c(2:9)],n_neighbors = 100, scale=T ,spread=1.5,n_epochs = 50)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
mycolors <- colorRampPalette(brewer.pal(9, "Set3"))(9)
p3 = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.8, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="bottom", byrow=F, reverse=TRUE)) + 
  xlab("umap1") + ylab("umap2") + theme(legend.position="none") + xlim(-10, 10) + ylim(-10, 10)
print(p3)

################### COLORED BY MARKER ############################################################################


for (i in c(2:9)){
  expression <- df[, i]
  p = ggplot(uplot,aes(x,y, color=expression))+ geom_point(size=0.8, stroke=0, alpha=0.7)+
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
  
}


##################################################################################################################
########################## TUMOR CELL UMAPS ######################################################################

#read.csv("TMA_annotated_single_cell_data.csv")

all_celltypes_24092020 <- all_celltypes_tumor
all_celltypes_24092020$Subtype <- as.character(all_celltypes_24092020$Subtype)

all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Hyperfunctional epithelial"), "Subtype"] <- "Functional epithelial"

#all_cell_types_13102020_tumor <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
#all_cell_types_13102020_tumor <- all_cell_types_13102020_tumor[-which(all_cell_types_13102020_tumor$GlobalCellType == "Negative"),]

df <- all_celltypes_24092020[, c("Sample","cCasp3", "pSTAT1", "Ki67",
                                 "PDL1","vimentin","Ecadherin","P21", "Subtype")]
df$Subtype <- as.factor(df$Subtype)


cell_type <- df$Subtype

umap_s = uwot::umap(df[,c(2:8)],n_neighbors = 80, scale=T ,spread=0.7,n_epochs = 30)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
mycolors <- rev(colorRampPalette(brewer.pal(7, "Accent"))(7))
p4 = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="bottom", byrow=F, reverse=TRUE)) + 
  xlab("umap1") + ylab("umap2") + theme(aspect.ratio=1, legend.position="none") 
print(p4)


############################################### COLORED BY EXPRESSION ######################################

for (i in c(2:8)){
  expression <- df[, i]
  p = ggplot(uplot,aes(x,y, color=expression))+ geom_point(size=0.7, stroke=0, alpha=0.7)+
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
}






