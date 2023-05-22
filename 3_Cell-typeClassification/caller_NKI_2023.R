#Cell type caller used for NKI project
library(corrplot)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
theme_set(theme_pubr())
library(mclust)

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
# Load the external gates.R file that contains cell type markers
source('gates_NKI_2023.R')
source('cell_type_caller_functions.R')

#File that describe association of cores with patients IDs
cycif2samples <- read.table(file="D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/3_Cell-typeClassification/tCicyf_patients-cores_pass.csv", sep=",", header = TRUE)

#List with the cell-types for subsequent cell.type calling
#A for loop will iterate through the this list, it will take the cells with the corresponding cell.labels then
#perform the Tribus-gating using the tribus.gate, and return the labels to the global_celltypes
celltypes.gating <- list()
celltypes.gating[["Lymphoid"]] <- list(
  cell.labels=c("Lymphoids1","Lymphoids2"),
  tribus.gate=Lymphoid.gate)
celltypes.gating[["CD8"]] <- list(
  cell.labels=c("CD8.T.cells"),
  tribus.gate=cd8.gates)
celltypes.gating[["Myeloid"]] <- list(
  cell.labels=c("Myeloids1","Myeloids2","Myeloids3"),
  tribus.gate=Myeloid.gate)

######################################################################################
#Reading information for samples (patient ids) that signal for CK7 and ECadherin needs to be gated
#This Gating threshold were manually selecting to improve the classification of cancer and immune cells
to.gate.cancer.cells <- read.table(file="D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/utils/Gates_ECadherin_CK7_patients.csv",
                                   sep=",", header = TRUE)



######################################################################################
#Reading cores to ignore
cores2ignore1 <- read.table(file="D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/utils/Total_cores_to_ignore.csv", sep=",")
cores2ignore2 <- read.table(file="D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/utils/Total_extra-cores_to_ignore.csv", sep = ",", header = TRUE)
cores2ignore3 <- read.table(file="D:/users/fperez/NKI_TMAs_AF/devNKI-scripts/utils/Total_cores_to_ignore_extra3.csv", sep = ",", header = TRUE)

cores2ignore3 <- cores2ignore3[,c(1,2)]
cores2ignore2 <- cores2ignore2[,-3]
colnames(cores2ignore1) <- colnames(cores2ignore2)
cores2ignore <- rbind(cores2ignore2, cores2ignore1, cores2ignore3)

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
outputfolders <- "Tribus_Cell-segment_202306/"
segmentation.suffix <- "_Probabilities_cytomask2" #This is to make it compatible with Napary
coresCoordsPath <- "dearray/cropCoords/"
coresCoordsSuffix <- "_cropCoords.csv"
Pixels.out.core <- 1175

# List all the data
files <- list.files(paste0(project.folder,input.folder), full.names = T)
files <- files[grepl(pattern = input.prefix, x=files)]
files <- files[grepl(pattern = input.suffix, x=files)]
set.seed(1) #For reproducibility 


for (f in files){
        print(paste0("Analyzing file: ", f))
        ########################Do this for every file #########################################
        print("Reading input data")
        data <- read.csv(f, sep=',', stringsAsFactors = F, header=TRUE)
        slide <- gsub(input.prefix, '',basename(f))
        slide <- gsub(input.suffix, '',slide)
        
        #List of cores by tumor sample (patient.id)
        cores.patient <- cycif2samples[cycif2samples$cycif.slide %in% slide,]
        cores.patient$cycif.core.id <- as.numeric(gsub("core","",cores.patient$cycif.core.id))
        
        #There are controls that should not be used for the normalization
        #Ignoring cores that should not be used for subsequent analysis
        cores2ignore.slide <- cores2ignore[cores2ignore$Slide %in% slide,]
        cores.patient <- cores.patient[!cores.patient$cycif.core.id %in% cores2ignore.slide$Core,]
        
        cores.to.ignore <- unique(data$CoreId[which(!data$CoreId %in% cores.patient$cycif.core.id)])
        if (length(cores.to.ignore) > 0){
           print("Next cores are ignored:")
           print(cores.to.ignore)
           data <- data[which(data$CoreId %in% cores.patient$cycif.core.id),]
        }
        
        outout.folder.name <- paste0(project.folder, "/", slide, "/", outputfolders)
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
        
        
        ############################### Plots to show summaries ######################################
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
        
        ###### Identification of cancer cells according to CK7, ECadherin, Vimentin and aSMA signal normalized by sample (patient) #############
        ##########Then further detection of cancer cells using a two Gaussian mixture approach for some samples
        ##########For some other samples with atypical signal, the signal was gated according to CK7, ECadherin, Vimentin and aSMA
        cores.to.gate <- cores.patient[cores.patient$patient %in% to.gate.cancer.cells$Patient,"cycif.core.id"]
        patients.to.gate <- cores.patient[cores.patient$patient %in% to.gate.cancer.cells$Patient,]
        patients.to.gaussian <- cores.patient[!cores.patient$patient %in% to.gate.cancer.cells$Patient,]
        data.filtered.to.gate <- data.filtered[data.filtered$CoreId %in% cores.to.gate,]
        data.filtered.to.gaussian <- data.filtered[!data.filtered$CoreId %in% cores.to.gate,]
        cancer.labels.gaussian <- gaussian.mix.cancer.cell.detection(data.filtered.to.gaussian, patientIDstable = patients.to.gaussian)
        cancer.labels.gate <- gating.cancer.cell.detection(data.filtered.to.gate, patientIDstable = patients.to.gate,
                                                               gating.file=to.gate.cancer.cells)
        is.Cancer.cell <- rep("NA", nrow(data.filtered)) 
        data.filtered.cancer.label <- cbind(data.filtered[,c(1:3)], is.Cancer.cell)
        data.filtered.cancer.label[data.filtered.cancer.label$CoreId %in% cores.to.gate, "is.Cancer.cell"] <- cancer.labels.gate$is.Cancer.cell
        data.filtered.cancer.label[!data.filtered.cancer.label$CoreId %in% cores.to.gate, "is.Cancer.cell"] <- cancer.labels.gaussian$is.Cancer.cell
        rm(is.Cancer.cell, cancer.labels.gate, cancer.labels.gaussian, data.filtered.to.gate, data.filtered.to.gaussian)
        
        
        #####For all the dataset trim channels to 0.999 & 0.001 percentile
        data.filtered.trim <- z.trimming(data.filtered[,cols])
        data.filtered[,cols] <- data.filtered.trim[[1]]
        
        ###Selecting non-cancer cells
        data.filtered.non.cancer <- data.filtered[which(data.filtered.cancer.label$is.Cancer.cell == "NO"),]
        
        #Detecting cells with an expression of all channels lower than the median channel intensity per slide
        #Ignoring those detected cells and saving those in vector any.high.expression
        #Selecting only those with expression of any channel above the median
        colnames.global <- unique(unlist(global.gates))
        colnames.global <- colnames.global[which(!colnames.global %in% 'Eccentricity')]
        data.filtered.non.cancer.scaled <- as.matrix(scale(data.filtered.non.cancer[,colnames.global]))
        quantile.by.channel <- sapply(colnames.global,function(x){
          quantile(data.filtered.non.cancer.scaled[,x],0.25)
        })
        any.high.expression <- sapply(1:nrow(data.filtered.non.cancer.scaled), function(x){
          any(data.filtered.non.cancer.scaled[x,] >= quantile.by.channel)
        })
        data.filtered.non.cancer.sel <- data.filtered.non.cancer[which(any.high.expression),]
        
        ############  Tribus, cell-type call ########################
        ####Cell-type calling for non-cancer cells ###########
        print("Running classifier, first gate")
        gating_results <- cellTypeCaller(data.filtered.non.cancer.sel, global.gates, "GlobalCellType", folder.name = outout.folder.name,
                                         grid.size.xdim=20, hierarchical.trees = FALSE, scaling = TRUE)
        
        #Merging results with signal data
        globalTypes <- gating_results[[1]] #Cell_types labels
        colnames(globalTypes) <- c("CellId2","GlobalCellType")
        
        
        #Cells with low intensity for all channels are labeled as "Other"
        GlobalCellTypes <- rep("Others",nrow(data.filtered.non.cancer))
        GlobalCellTypes[which(any.high.expression)] <- globalTypes$GlobalCellType
        
        #Gatting B-cells
        GlobalCellTypes[ which(as.vector(scale(data.filtered.non.cancer$CD20)) >= 3)] <- "B.cells"
        
        #Merging results with previous cancer-cell assignment
        data.filtered.types <- data.filtered
        data.filtered.types$GlobalCellType <- NA
        data.filtered.types$GlobalCellType[which(data.filtered.cancer.label$is.Cancer.cell == "YES")] <- "Cancer"
        data.filtered.types$GlobalCellType[which(data.filtered.cancer.label$is.Cancer.cell == "NO")] = GlobalCellTypes
        
        print("Cell types proportions found at global gate:")
        print((table(data.filtered.types$GlobalCellType) / nrow(data.filtered.types)) * 100)
        
        ##########Second call, immune cell types

        ###########Performing subsequent Immune cell type.calling
        #Scaling immune markers just for the detected immune cells; but for cancer markers the scaling was done using all cells
        cells.for.gate <- data.filtered.types[which(!data.filtered.types$GlobalCellType %in% c("Cancer","Stromal1","Stromal2","Stromal3","Stromal4", "Others","B.cells")), ]
        cells.for.gate.data  <- cells.for.gate[,colnames(cells.for.gate) %in% unique(unlist(immune.gates))]

        #cellTypeCaller
        Immune.Types <- cellTypeCaller(cells.for.gate, immune.gates, "ImmuneCelltype", folder.name = outout.folder.name, hierarchical.trees=FALSE, scaling = TRUE)[[1]]
        Immune.Types <- Immune.Types[,2]
        Immune.Types2   <- cbind(cells.for.gate, Immune.Types)
        numbers_celltypes <- table(Immune.Types2$Immune.Types)
        percentages_celltypes <- table(Immune.Types2$Immune.Types) * 100 / nrow(Immune.Types2)
        print(percentages_celltypes)

        #Calling sub-types of immune cells
        input.subgating <- Immune.Types2
        for (j in 1:length(celltypes.gating)){
          cells.subgate <- input.subgating[input.subgating$Immune.Types %in% celltypes.gating[[j]]$cell.labels,]
          subgate.types <- cellTypeCaller(cells.subgate, celltypes.gating[[j]]$tribus.gate, names(celltypes.gating)[j],
                                          folder.name = outout.folder.name, grid.size.xdim=5, hierarchical.trees=FALSE)[[1]]
          #Merging results in corresponding column
          input.subgating[input.subgating$Immune.Types %in% celltypes.gating[[j]]$cell.labels,"Immune.Types"] <- subgate.types[,2]
        }
        Immune.Types3 = input.subgating
        #Cells with high proportion of CD15 are called CD15.MY
        Immune.Types3[ which(as.vector(scale(cells.for.gate$CD15)) >= 3),"Immune.Types"] <- "CD15.MY"
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
        umaps.plots <- channel.UMAPs(immune.data, cell_type, sub.sample = TRUE, n.sampling = 10000)
        ggsave(file=paste0(outout.folder.name,"Immune_UMAPs_channels_untruncated.pdf"), arrangeGrob(grobs = umaps.plots[[1]], ncol = 4), height = 14, width = 16)
        ggsave(file=paste0(outout.folder.name,"Immune_UMAPs_channels_untruncated.png"), arrangeGrob(grobs = umaps.plots[[1]], ncol = 4), height = 5.3, width = 6.1)
        ggsave(file=paste0(outout.folder.name,"Immune_UMAPs_celltypes_untruncated.pdf"), umaps.plots[[2]], height = 6.89, width = 7.93)
        detach("package:ComplexHeatmap", unload = TRUE) #This avoid conflict with the cellTypeCaller function


        ##########Concatenating the results from immune cells to the initial table of global.gates
        df_merged_global.aux <- df_merged_global[,-ncol(df_merged_global)] #Ignoring last column
        remaining.cells <- data.filtered.types[which(data.filtered.types$GlobalCellType %in% c("Cancer","Stromal1","Stromal2","Stromal3","Stromal4","Others","B.cells")),]
        data.filtered.types2 <- rbind(df_merged_global.aux, remaining.cells)

        ####Merge with initial cell-table, and add labels for those cells ignored before the classification
        selected.cols <- c("Row_number", "GlobalCellType") #Row_number is the the first column's name
        data.filtered.types.cols  <- data.filtered.types2[,selected.cols]
        data.raw.types <- merge(data, data.filtered.types.cols, all.x = TRUE, by="Row_number")
        data.raw.types$GlobalCellType[which(is.na(data.raw.types$GlobalCellType))] <- "SignalQC"
        data.raw.types$GlobalCellType[dapi.low.chanel.cells] <- "LowDapi"
        data.raw.types$GlobalCellType[which(data.raw.types$lost == "True")] <- "LostCell"
        data.raw.types$GlobalCellType[which(data.raw.types$Inside_ROI == 1)] <- "InsideROI"
        data.raw.types$out_of_core <- NA #Label if the cell correspond to the current core

        percentages_celltypes <- round(table(data.raw.types$GlobalCellType) * 100 / nrow(data.raw.types),4)
        print(paste0("Final cell types proportions for ", slide, ":"))
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

        # ###Code for saving and for checking the cores outside the bounding box
        print(paste0("Saving output files in folder: ", outout.folder.name))
        for (core.id in unique(data.raw.types$CoreId)){
                core <- data.raw.types[data.raw.types$CoreId == core.id,c(2:ncol(data.raw.types))]
                names(core)[which(names(core) == "CellId")] <- "Cellid"
                names(core)[which(names(core) == "CoreId")] <- "Core_Names"
                #Calculating the distance of each cell to the center of the core
                cropbox = read.table(file=paste0(project.folder, slide, "/", coresCoordsPath, core.id, coresCoordsSuffix), sep=",")
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

        write.table(data.raw.types, file=paste0(outout.folder.name, slide, "_cellTypes.csv"), sep=",", row.names=FALSE)
}
