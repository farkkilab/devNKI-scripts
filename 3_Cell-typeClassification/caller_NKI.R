#Cell type caller used for NKI project
library(corrplot)
library(reshape2)
library(ggplot2)

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


#Load functions
setwd("D:/users/fperez/NKI_TMAs_AF/R_tribus/")
source('2b_cell_type_caller.R')
source('qc_funcs.R')
# Load the external gates.R file that contains cell type gatings
source('gates_NKI.R')


############################################################################################
########## Z-score outlier detection and normalization function ############################
############################################################################################

#As input a data.frame with only intensity signal, transform it to log2
#It will normalize the data by z-score
#Then will remove values with z-score above 6
#For values with z-score above 4, it will give the same values as the intensity of 4
#It will also return the columns removed from your original dataset

z.trimming <- function(signal.data, z.cutoff = 7, z.for.max.intensity = 4){
        z.scores.data <- lapply(1:ncol(signal.data), function(x){
                z.scores <- (signal.data[,x] - mean(signal.data[,x])) / sd(signal.data[,x])
                return(z.scores)
        })
        #Adjusting the data frame
        to.remove.pos.all <- NULL
        for (i in (1:length(z.scores.data))){
                #First check if there are outliers
                if (length(which(z.scores.data[[i]] >= z.for.max.intensity)) > 0) {
                        max.value <- z.for.max.intensity * sd(signal.data[[i]]) + mean(signal.data[[i]])
                        #to.adjust.pos <- which(z.scores.data[[i]] >= z.for.max.intensity & z.scores.data[[i]] < z.cutoff)
                        to.adjust.pos <- which(z.scores.data[[i]] >= z.for.max.intensity)
                        signal.data[to.adjust.pos,i] <- max.value
                        #to.remove.pos <- which(z.scores.data[[i]] >= z.cutoff)
                        #to.remove.pos.all <- c(to.remove.pos.all, to.remove.pos)
                }
        }
        #Remove those values selected
        if (length(to.remove.pos.all) > 0){
                signal.data <- signal.data[-to.remove.pos.all,]
        }
        result <- list(signal.data, to.remove.pos.all)
        return(result)
}
############################################################################################
############################################################################################

# Edit these to correspond to the naming of your dna and background channels
# E.g my DNA channels were Hoechst_1 etc. so a regexp to catch that is Hoech
# Separate bg channels with a pipe. Here's a couple extras in case you need them ||||||||||||||||||||||||||
#dnachannel = 'Dapi'
#bgchannels = 'Background1|Background2|Background3'

## Set the directory to your working directory that contains the output data from QC
#dir.create("plots", showWarnings = FALSE)
#dir.create("gated", showWarnings = FALSE)

###Prefix and suffix for input cell quantification tables after QC
input.prefix <- "annotated_"
input.suffix <- ".csv"
#I ussed the next file name structure "annotated_TMA_18_810.csv",
#all before suffix and after prefix is the sample name

###Input and output folders
project.folder <- "D:/users/fperez/NKI_TMAs_AF/"
input.folder <- "Cell_QCs4/"
outputfolders <- "Tribus_celltype2/"

# List all the data
files <- list.files(paste0(project.folder,input.folder), full.names = T)
files <- files[grepl(pattern = input.prefix, x=files)]

for (f in files){
 f <- files[1]
        print(paste0("Analyzing file: ", f))
        ########################Do this for every file #########################################
        print("Reading input data")
        data <- read.csv(f, sep=',', stringsAsFactors = F, header=TRUE)
        sample <- gsub(input.prefix, '',basename(f))
        sample <- gsub(input.suffix, '',sample)

        ############################## Some QCs before gating ###################################
        ###Log2 transformation of signal intensity
        print("Performing data normalization")
        cols <- which(colnames(data) =='DNA1'):c(which(colnames(data) =='Area') - 1)
        data.filtered = data
        data.filtered[,cols] = log2(data[,cols])
        
        
        
        ###Remove cells with very low expression of DAPI
        dapi.low.chanel.cells <- which(data.filtered$DNA1 < 10.5)
        print(paste0("Cells with very low dapi: ", length(dapi.low.chanel.cells)))
        data.filtered <- data.filtered[-dapi.low.chanel.cells,]
        
        #Ignore cells labeled as lost by the DAPI concordance between channels
        data.filtered <- data.filtered[(data.filtered$lost == "False"),]
        
        ##Trim high values, and remove extremely high values
        gate.cols <- which(colnames(data.filtered) %in% unique(unlist(global.gates)))
        data.filtered.trim <- z.trimming(data.filtered[,gate.cols])
        #data.filtered.ad <- data.filtered[-data.filtered.trim[[2]],]
        #print(paste0("Cells with extremely high signal: ", length(data.filtered.trim[[2]])))
        #data.filtered.ad[,gate.cols] <- data.filtered.trim[[1]]
        data.filtered[,gate.cols] <- data.filtered.trim[[1]]
        
        #Remove the 1% cells with signal so high in the channels used for global gating
        #channels.interest <- unique(unlist(global.gates))
        #high.cells  <- lapply(channels.interest, function(x){
        #  cutvalue <- quantile(data.filtered[,x], probs = c(0.999))
        #  which(data.filtered[,x] >= cutvalue)
        #})
        #high.cells.to.remove <- unique(unlist(high.cells))
        #print(paste0("Cells with very high global gates channel intensity: ", length(high.cells.to.remove)))
        #data.filtered <- data.filtered[-c(high.cells.to.remove),]
        
        #Check density distribution for of some channels
        #plot(density(data.filtered$CD20))
        #plot(density(data.filtered$DNA1))
        
        ################################  Gate calls #############################################
        
        #Output folder
        outout.folder.name <- paste0(project.folder, "/", sample, "/", outputfolders)
        dir.create(outout.folder.name)
        
        # First call global cell types
        print("Running classifier")
        globalTypes <- cellTypeCaller(data.filtered, global.gates, "GlobalCellType", folder.name = outout.folder.name)
        
        colnames(globalTypes) <- c("CellId2","GlobalCellType")
        data.filtered.types   <- cbind(data.filtered, globalTypes)
        selected.cols <- c("X", "GlobalCellType")
        data.filtered.types.cols  <- data.filtered.types[,selected.cols]
        
        data.raw.types <- merge(data, data.filtered.types.cols, all.x = TRUE, by="X")
        data.raw.types$GlobalCellType[which(is.na(data.raw.types$GlobalCellType))] <- "SignalQC"
        data.raw.types$GlobalCellType[dapi.low.chanel.cells] <- "LowDapi"
        data.raw.types$GlobalCellType[which(data.raw.types$lost == "True")] <- "LostCell"
        
        print("Cell types proportions found:")
        print((table(data.raw.types$GlobalCellType) / nrow(data.raw.types)) * 100)
        
        print(paste0("Saving output files in folder:", outout.folder.name))
        for (core.id in unique(data.raw.types$CoreId)){
                core <- data.raw.types[data.raw.types$CoreId == core.id,c(2:ncol(data.raw.types))]
                names(core)[which(names(core) == "CellId")] <- "Cellid"
                names(core)[which(names(core) == "CoreId")] <- "Core_Names"
                core.set <- core[,c("Core_Names","Cellid","lost","GlobalCellType")]
                core.set[,1] <- rep(paste0("core",core.id, "_Probabilities_cytomask2"), nrow(core.set))
                write.table(core.set, file=paste0(outout.folder.name,"core", core.id,"_cellTypesNapari.csv"),sep=",", row.names=FALSE)
        }
        write.table(data.raw.types, file=paste0(outout.folder.name, sample, "core", "_cellTypes.csv"), sep=",", row.names=FALSE)
}


