#Cell type caller used for NKI project
library(corrplot)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
theme_set(theme_pubr())

# For these you need Biocmanager if not installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("FlowSOM")
library(flowCore)
library(FlowSOM)

library(sjstats)
library(BBmisc)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(foreach)
library(doParallel)
library(uwot)

library(future.apply)
library(doFuture)
detach("package:ComplexHeatmap", unload = TRUE)

heatmaps <- list()


#Load functions
setwd("D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/3_Cell-typeClassification/")
source('2b_cell_type_caller.R')
# Load the external gates.R file that contains cell type gatings
source('gates_NKI.R')
source('cell_type_caller_functions.R')


#List with the cell-types for subsequent gatings
#A for loop will iterate throug the this list, take the cells with the corresponding cell.labels then
#perform the Tribus-gating using the tribus.gate, and return the labels to the global_celltypes
celltypes.gating <- list()
celltypes.gating[["Lymphoid"]] <- list(
  cell.labels=c("Lymphoids1","Lymphoids2"),
  tribus.gate=Lymphoid.gate)
celltypes.gating[["Myeloid"]] <- list(
  cell.labels=c("Myeloids1","Myeloids2","Myeloids3"),
  tribus.gate=Myeloid.gate)



####################################################################################################
####################################################################################################

#I used the next file name structure "annotated_TMA_18_810_ROIs-label.csv"
#Last column of csv files have information binary information (0 or 1). Value of 1 if the sample was in a ROI for cell-removal
#Second last column have information if a cell was lost using the DAPI thresholds in cycif-suit
#TMA_18_810 correspond to the slide name

###Prefix and suffix for input cell quantification tables after QC
input.prefix <- "annotated_"
input.suffix <- "_ROIs-label.csv"


###Input and output folders
project.folder <- "D:/users/fperez/NKI_TMAs_AF/"
input.folder <- "Cell-segment2_QC/"
outputfolders <- "Tribus_Cell-segment2/"
segmentation.suffix <- "_Probabilities_cytomask2" #This is to make it compatible with Napary
coresCoordsPath <- "dearray/cropCoords/"
coresCoordsSuffix <- "_cropCoords.csv"
Pixels.out.core <- 1200

# List all the data
files <- list.files(paste0(project.folder,input.folder), full.names = T)
files <- files[grepl(pattern = input.prefix, x=files)]
files <- files[grepl(pattern = input.suffix, x=files)]


for (f in files){
        print(paste0("Analyzing file: ", f))
        ########################Do this for every file #########################################
        print("Reading input data")
        data <- read.csv(f, sep=',', stringsAsFactors = F, header=TRUE)
        sample <- gsub(input.prefix, '',basename(f))
        sample <- gsub(input.suffix, '',sample)

        outout.folder.name <- paste0(project.folder, "/", sample, "/", outputfolders)
        dir.create(outout.folder.name)

        ############################## Some QCs before gating ###################################
        ###Log2 transformation of signal intensity
        print("Performing data normalization")
        cols <- which(colnames(data) =='DNA1'):c(which(colnames(data) =='Area') - 1)
        data.filtered = data
        data.filtered[,cols] = log2(data[,cols])


        #####Ignore cells with very low expression of DAPI
        dapi.low.chanel.cells <- which(data.filtered$DNA1 < 10.5)
        print(paste0("Cells with very low dapi: ", length(dapi.low.chanel.cells)))
        data.filtered <- data.filtered[-dapi.low.chanel.cells,]

        #####Ignore cells labeled as lost by the DAPI concordance between channels and in ROIs
        data.filtered <- data.filtered[(data.filtered$lost == "False"),]
        data.filtered <- data.filtered[(data.filtered$Inside_ROI == "0"),]

        #####Trim channels to 0.999 & 0.001 percentile
        data.filtered.trim <- z.trimming(data.filtered[,cols])
        data.filtered[,cols] <- data.filtered.trim[[1]]

        ##################################  Gate calls ###############################################

        ###########Signal intensity by channel for the global gates
        source('qc_functions.R')
        densities.by.channels <- channel.densities(data.filtered[,cols])
        ggsave(file=paste0(outout.folder.name,"GlobalCellType_channel_densities.pdf"), arrangeGrob(grobs = densities.by.channels, ncol = 4), height = 25, width = 16)
        ggsave(file=paste0(outout.folder.name,"GlobalCellType_channel_densities.png"), arrangeGrob(grobs = densities.by.channels, ncol = 4), height = 6.89, width = 7.93)
        detach("package:ComplexHeatmap", unload = TRUE) #This avoid conflict with the cellTypeCaller function


        colnames.interest <- colnames(data.filtered[,cols])[c(!(grepl('DNA', colnames(data.filtered[,cols]))
                                                               | grepl('BG', colnames(data.filtered[,cols]))))]
        pdf(paste0(outout.folder.name,"GlobalCellType_channel_correlations_all.pdf"))
        heatmap(cor(data.filtered[,colnames.interest]))
        dev.off()

        not.celltype.channels <- c("pTBK1","PDL1_488","CD45RO","pSTAT1", "PDL1_2",
                                   "PDL1", "PDL1_555", "PD1", "FOXP3", "LaminB1",
                                   "MHCI",  "yH2AX", "cPARP1", "GranzymeB", "TIM3", "Ki67")

        colnames.interest2 <- colnames.interest[!colnames.interest %in% not.celltype.channels]
        pdf(paste0(outout.folder.name,"GlobalCellType_channel_correlations_markers.pdf"))
        heatmap(cor(data.filtered[,colnames.interest2]))
        dev.off()

        ##########First call global cell types
        print("Running classifier, first gate")
        source('gates_NKI.R')
        #gating_results <- cellTypeCaller(data.filtered, global.gates, "GlobalCellType", folder.name = outout.folder.name, grid.size.xdim=25, hierarchical.trees = FALSE)
        gating_results <- cellTypeCaller(data.filtered, global.gates, "GlobalCellType", folder.name = outout.folder.name, grid.size.xdim=20, hierarchical.trees = FALSE)

        globalTypes <- gating_results[[1]] #Cell_types labels
        nodes.scores <- gating_results[[2]] #To get the scores of the nodes, and find the cells with higher cancer score
        top.cancer.nodes <- order(nodes.scores[,"Cancer"], decreasing = TRUE)[1:2]
        nodes.by.cell <- gating_results[[3]] #The node in which corresponds each cell

        #Merging results with signal data
        colnames(globalTypes) <- c("CellId2","GlobalCellType")
        data.filtered.types   <- cbind(data.filtered, globalTypes)

        print("Cell types proportions found at global gate:")
        print((table(data.filtered.types$GlobalCellType) / nrow(data.filtered.types)) * 100)

        source('qc_functions.R')
        ind.plots = index_plots(data.filtered.types,data.filtered.types, gates.input=c(global.gates), scaling = "z.score")
        pdf(file= paste0(outout.folder.name,"Global-gate_celltypes_markerExpression.pdf"), height = 8, width = 10)
        print(ind.plots[[3]])
        dev.off()
        detach("package:ComplexHeatmap", unload = TRUE) #This avoid conflict with the cellTypeCaller function

        ##########Third call, immune cell types

        ###########Performing subsequent Immune gating
        #Scaling immune markers just for the detected immune cells; but for cancer markers the scaling was done using all cells
        cells.for.gate <- data.filtered.types[which(!data.filtered.types$GlobalCellType %in% c("Cancer","Stromal1","Stromal2", "Others")), ]
        cells.for.gate.data  <- cells.for.gate[,colnames(cells.for.gate) %in% unique(unlist(immune.gates))]

        pdf(paste0(outout.folder.name,"ImmuneCellType_channel_correlations2.pdf"))
        heatmap(cor(cells.for.gate.data))
        dev.off()

        ######Signal intensity by channel
        densities.by.channels <- channel.densities(cells.for.gate.data)
        ggsave(file=paste0(outout.folder.name,"Immune_channel_densities.pdf"), arrangeGrob(grobs = densities.by.channels, ncol = 4), height = 14, width = 16)
        ggsave(file=paste0(outout.folder.name,"Immune_channel_densities.png"), arrangeGrob(grobs = densities.by.channels, ncol = 4), height = 6.89, width = 7.93)

        #cellTypeCaller
        Immune.Types <- cellTypeCaller(cells.for.gate, immune.gates, "ImmuneCelltype", folder.name = outout.folder.name, hierarchical.trees=FALSE, scaling = TRUE)[[1]]
        Immune.Types <- Immune.Types[,2]
        Immune.Types2   <- cbind(cells.for.gate, Immune.Types)
        numbers_celltypes <- table(Immune.Types2$Immune.Types)
        percentages_celltypes <- table(Immune.Types2$Immune.Types) * 100 / nrow(Immune.Types2)
        print(percentages_celltypes)

        write.table(t(numbers_celltypes), file=paste0(outout.folder.name,"ImmuneCellType_cells-numbers.csv"), sep=",", row.names = FALSE)
        write.table(t(percentages_celltypes), file=paste0(outout.folder.name,"ImmuneCellType_cells-percentage.csv"), sep=",", row.names = FALSE)

        input.subgating <- Immune.Types2
        for (j in 1:length(celltypes.gating)){
          cells.subgate <- input.subgating[input.subgating$Immune.Types %in% celltypes.gating[[j]]$cell.labels,]
          subgate.types <- cellTypeCaller(cells.subgate, celltypes.gating[[j]]$tribus.gate, names(celltypes.gating)[j],
                                          folder.name = outout.folder.name, grid.size.xdim=5, hierarchical.trees=FALSE)[[1]]
          #Merging results in corresponding column
          input.subgating[input.subgating$Immune.Types %in% celltypes.gating[[j]]$cell.labels,"Immune.Types"] <- subgate.types[,2]
        }
        Immune.Types3 = input.subgating
        Immune.Types3[Immune.Types3$Immune.Types == "CancerImmune","Immune.Types"] <- "Cancer"

        table(Immune.Types3$Immune.Types) * 100 /nrow(Immune.Types2)

        ######Performing QC plots for the immune gating
        source('qc_functions.R')
        GlobalCellType = "GlobalCellType" #To change variables inside the qc_functions.R
        cell.type.column = "GlobalCellType"  #To change variables inside the qc_functions.R
        ImmuneCellType = "Immune.Types"   #To change variables inside the qc_functions.R
        #Adding the immune gate labels to the global gate column
        df_merged_global = combine_gates(Immune.Types3, immune_labels = unique(Immune.Types3$GlobalCellType))
        #pie_charts_abundances  = all_pie_charts(df_merged_global, group.by = "CoreId")

        #Immune.plots.data <- df_merged_global[df_merged_global$GlobalCellType != "Cancer",]

        ind.plots = index_plots(df_merged_global,df_merged_global, gates.input=c(immune.gates,Lymphoid.gate, Myeloid.gate), scaling = "z.score")
        pdf(file= paste0(outout.folder.name,"Immune_celltypes_markerExpression.pdf"))
        print(ind.plots[[3]])
        dev.off()

        ##########Generating UMAPs for immune
        immune.interesting.channels <- unique(unlist(c(Myeloid.gate, immune.gates, Lymphoid.gate)))
        immune.data <- df_merged_global[df_merged_global$GlobalCellType != "Cancer" ,immune.interesting.channels]
        immune.types <- df_merged_global[df_merged_global$GlobalCellType != "Cancer" ,"GlobalCellType"]

        #cell_type <- factor(Immune.Types, levels=(unique(Immune.Types))) #Cell type for each row
        cell_type <- factor(immune.types, levels=(unique(immune.types))) #Cell type for each row
        umaps.plots <- channel.UMAPs(immune.data, cell_type, sub.sample = TRUE, n.sampling = 20000)
        ggsave(file=paste0(outout.folder.name,"Immune_UMAPs_channels_untruncated.pdf"), arrangeGrob(grobs = umaps.plots[[1]], ncol = 4), height = 14, width = 16)
        ggsave(file=paste0(outout.folder.name,"Immune_UMAPs_channels_untruncated.png"), arrangeGrob(grobs = umaps.plots[[1]], ncol = 4), height = 5.3, width = 6.1)
        ggsave(file=paste0(outout.folder.name,"Immune_UMAPs_celltypes_untruncated.pdf"), umaps.plots[[2]], height = 6.89, width = 7.93)
        detach("package:ComplexHeatmap", unload = TRUE) #This avoid conflict with the cellTypeCaller function


        ##########Concatenating the results from immune cells to the initial table of global.gates
        df_merged_global.aux <- df_merged_global[,-ncol(df_merged_global)] #Ignoring last column
        remaining.cells <- data.filtered.types[which(data.filtered.types$GlobalCellType %in% c("Cancer","Stromal1","Stromal2","Others")),]
        data.filtered.types2 <- rbind(df_merged_global.aux, remaining.cells)

        ####Merge with initial cell-table, and add labels for those cells ignored before the classification
        selected.cols <- c("Row_number", "GlobalCellType") #Row_number is the the first column's name
        data.filtered.types.cols  <- data.filtered.types2[,selected.cols]
        data.raw.types <- merge(data, data.filtered.types.cols, all.x = TRUE, by="Row_number")
        data.raw.types$GlobalCellType[which(is.na(data.raw.types$GlobalCellType))] <- "SignalQC"
        data.raw.types$GlobalCellType[dapi.low.chanel.cells] <- "LowDapi"
        data.raw.types$GlobalCellType[which(data.raw.types$lost == "True")] <- "LostCell"
        data.raw.types$GlobalCellType[which(data.raw.types$Inside_ROI == 1)] <- "InsideROI"
        data.raw.types$GlobalCellType[which(data.raw.types$Row_number %in% others.ids)] <- "Others"
        data.raw.types$out_of_core <- NA #Label if the cell correspond to the current core

        plot(density(log2(data.raw.types[data.raw.types$GlobalCellType == "Others","CD8a"])))
        plot(density(log2(data.raw.types[data.raw.types$GlobalCellType == "CD8.T.cells","CD8a"])))

        percentages_celltypes <- round(table(data.raw.types$GlobalCellType) * 100 / nrow(data.raw.types),4)
        print(paste0("Final cell types proportions for ", sample, ":"))
        print(percentages_celltypes)

        data.raw.types.plots <- data.raw.types[!data.raw.types$GlobalCellType %in%  c("LostCell", "Others", "LowDapi", "InsideROI"),]

        source('qc_functions.R')
        pdf(file= paste0(outout.folder.name,"Total_celltypes_markerExpression_others.pdf"), height = 8, width = 10)
        ind.plots = index_plots(data.raw.types,data.raw.types, gates.input=c(global.gates, immune.gates, Lymphoid.gate, Myeloid.gate), scaling = "z.score", log2.transform = FALSE)
        print(ind.plots[[3]])
        dev.off()
        pdf(file= paste0(outout.folder.name,"Total-nonQC_celltypes_markerExpression.pdf"), height = 8, width = 10)
        ind.plots = index_plots(data.raw.types.plots,data.raw.types.plots, gates.input=c(global.gates[(names(global.gates) != "Others")], immune.gates, Lymphoid.gate, Myeloid.gate))
        print(ind.plots[[1]])
        print(ind.plots[[3]])
        dev.off()
        pdf(file= paste0(outout.folder.name,"Total-nonQC_celltypes_markerExpression_min-max.pdf"), height = 8, width = 10)
        ind.plots = index_plots(data.raw.types.plots,data.raw.types.plots, gates.input=c(global.gates[(names(global.gates) != "Others")], immune.gates, Lymphoid.gate, Myeloid.gate), scaling = "min.max")
        print(ind.plots[[1]])
        print(ind.plots[[3]])
        dev.off()
        detach("package:ComplexHeatmap", unload = TRUE) #This avoid conflict with the cellTypeCaller function

        ###Code for saving and for checking the cores outside the bounding box
        print(paste0("Saving output files in folder: ", outout.folder.name))
        for (core.id in unique(data.raw.types$CoreId)){
                core <- data.raw.types[data.raw.types$CoreId == core.id,c(2:ncol(data.raw.types))]
                names(core)[which(names(core) == "CellId")] <- "Cellid"
                names(core)[which(names(core) == "CoreId")] <- "Core_Names"
                #Calculating the distance of each cell to the center of the core
                cropbox = read.table(file=paste0(project.folder, sample, "/", coresCoordsPath, core.id, coresCoordsSuffix), sep=",")
                center.cropBox = data.frame(x=(cropbox$V3 - cropbox$V1)/2,
                                            y=(cropbox$V4 - cropbox$V2)/2)
                distances = sqrt((core$X_position - center.cropBox$x)^2 + (core$Y_position - center.cropBox$y)^2)
                core$out_of_core <-c("TRUE","FALSE")[match(distances >= Pixels.out.core, c('TRUE', 'FALSE'))]
                core[which(core$out_of_core == "TRUE"),"GlobalCellType"] <- "Out.of.core"
                data.raw.types[data.raw.types$CoreId == core.id,"out_of_core"] <- core$out_of_core
                #Taking subset data and writing in Napary format
                core.set <- core[,c("Core_Names","Cellid","lost","GlobalCellType")]
                core.set[,1] <- rep(paste0("core",core.id, segmentation.suffix), nrow(core.set))
                write.table(core.set, file=paste0(outout.folder.name,"core", core.id,"_cellTypesNapari_others.csv"),sep=",", row.names=FALSE)
        }
        #Excluding cells out of the cores, InsideROIs, lost cells and LowDapi
        data.raw.types <- data.raw.types[data.raw.types$out_of_core == "FALSE",]
        data.raw.types <- data.raw.types[data.raw.types$GlobalCellType != "InsideROI",]
        data.raw.types <- data.raw.types[data.raw.types$GlobalCellType != "LostCell",]

        write.table(data.raw.types, file=paste0(outout.folder.name, sample, "_cellTypes.csv"), sep=",", row.names=FALSE)
}


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

#Probably useful code for future development



# L <- data.filtered.types
# L.cancer <- L[L$GlobalCellType == "Cancer",]
#
#
# png(filename = paste0("CK7-ECadherin_",sample,"-immune-cells.png"))
# par(mfrow = c(1, 2))
# plot(density(cells.for.gate$CK7), main="CK7 - Immune cells")
# plot(density(cells.for.gate$ECadherin), main="ECadherin - Immune cells")
# dev.off()
#
# png(filename = paste0("CK7-ECadherin_",sample,"-immune-cells.png"))
# par(mfrow = c(1, 2))
# plot(density(L[L$GlobalCellType == "Cancer","CK7"]), main="CK7 - Cancer cells")
# plot(density(L[L$GlobalCellType == "Cancer","ECadherin"]), main="ECadherin - Cancer cells")
# dev.off()
#
#
#
#
# ggsave(fig, filename=paste0("Density_CK7-Ecadherin_",sample, "-", names.data.set[i], ".png"), width = 18, height = 12, units = "cm")
#
# cols.selected <- intersect(intersect(intersect(intersect(cols.1,cols.2),cols.3), cols.4),cols.5)
# cells.for.gate.data <- cells.for.gate.data[,cols.selected]
# D <- scale(as.matrix(cells.for.gate.data))
# cor.D <- cor(D)
# heatmap(cor.D)
# corrplot(cor.D, order="hclust")
# prc = prcomp(t(D))
# plot(prc$x[,1], prc$x[,2], col = "white", xlab = "PC1", ylab = "PC2")
# text(prc$x[,1], prc$x[,2], labels = colnames(D))
# plot(prc$x[,2], prc$x[,3], col = "white", xlab = "PC2", ylab = "PC3")
# text(prc$x[,2], prc$x[,3], labels = colnames(D))


# #Find cancer cells if those are from the same node
# cells.node.channels <- cbind(data.filtered.types, nodes.by.cell)
# rows.cancer <- Immune.Types2[which(Immune.Types2$Immune.Types == "Cancer"),"Row_number"]
# cancer.inside.immune <- cells.node.channels[which(cells.node.channels$Row_number %in% rows.cancer),]
