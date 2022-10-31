library(FlowSOM)
library(ConsensusClusterPlus)

#source("/mnt/d//users/fperez/NKI_TMAs_AF/devNKI-scripts/3_Cell-typeClassification/qc_functions.R")
source("/mnt/d//users/fperez/NKI_TMAs_AF/devNKI-scripts/3_Cell-typeClassification/cell_type_caller_functions.R")
source("/mnt/d//users/fperez/Programs/ConsensusClusterPlus/R/ConsensusClusterPlus.R")


table1 <- read.table(file="/mnt/d/users/fperez/NKI_TMAs_AF/Tables/Stromal_cells.csv", header=TRUE, sep=",")


#interesting.channels <- c("Ki67","Vimentin","pSTAT1", "MHCI", "ECadherin", "CK7",
#                          "Area","Eccentricity","Perimeter","Roundness")

interesting.channels <- c("pTBK1","Ki67","Vimentin", "pSTAT1",
                           "MHCI", "ECadherin", "aSMA", "CD31", 
                           "Area","Eccentricity","Perimeter","Roundness")


table.data  <- log2(table1[,interesting.channels[1:8]])
data.filtered.trim <- z.trimming(table.data)
table.data <- data.filtered.trim[[1]]
table.data <- cbind(table.data, table1[,interesting.channels[9:length(interesting.channels)]])

table.data.scaled <- scale(table.data)

n = 12000
to.sample <- sort(sample(c(1:nrow(table.data.scaled)), n))
table.data.scaled.set <- table.data.scaled[to.sample,]


fsomclust = function(data, k){
  fSOM <- FlowSOM(data, colsToUse = colnames(data), 
                  xdim = 25, ydim = 25, nClus = k)
  cluster.ids <- fSOM$map$mapping[,1]
  metacluster.ids <- fSOM$metaclustering[cluster.ids]
  print("Ending FlowSOM")
  return(metacluster.ids)
}

title="/mnt/d/users/fperez/NKI_TMAs_AF/Cell_class_analysis/Consensus_clustering_FlowSom-Stromal_Linux2/"

t.data.set <- t(table.data.scaled.set)

colnames(t.data.set) <- paste0("S",1:ncol(t.data.set))
rm(table1)

results = ConsensusClusterPlus2(t.data.set, clusterAlg="km",
                                maxK=12,reps=200, pItem=0.8,pFeature=1, title=title,
                                distance="euclidean", plot="pngBMP", verbose=TRUE,
                                writeTable=TRUE)
