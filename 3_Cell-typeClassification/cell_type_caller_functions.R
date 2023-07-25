###A function to truncate the dataset, using the 99.9 percentile
#Other commented lines in the script were added to make the truncation using percentiles and ignores cells with a Z score so high
z.trimming <- function(signal.data, z.cutoff = 10, z.for.max.intensity = 2,max.quantile=0.999, min.quantile=0.001){
  z.scores.data <- lapply(1:ncol(signal.data), function(x){
    z.scores <- (signal.data[,x] - mean(signal.data[,x])) / sd(signal.data[,x])
    return(z.scores)
  })
  #Adjusting the data frame
  to.remove.pos.all <- NULL
  for (i in (1:length(z.scores.data))){
    #First check if there are outliers
    if (length(which(z.scores.data[[i]] >= z.for.max.intensity)) > 0) {
      #max.value <- z.for.max.intensity * sd(signal.data[[i]]) + mean(signal.data[[i]])
      #min.value <- z.for.max.intensity * sd(signal.data[[i]]) - mean(signal.data[[i]])
      max.value <- quantile(signal.data[,i], max.quantile)
      min.value <- quantile(signal.data[,i], min.quantile)
      #to.adjust.pos <- which(z.scores.data[[i]] >= z.for.max.intensity & z.scores.data[[i]] < z.cutoff)
      #to.adjust.pos <- which(z.scores.data[[i]] >= z.for.max.intensity)
      #signal.data[to.adjust.pos,i] <- max.value
      #to.remove.pos <- which(z.scores.data[[i]] >= z.cutoff)
      #to.remove.pos.all <- c(to.remove.pos.all, to.remove.pos)
      to.adjust.pos <- which(signal.data[,i] >= max.value)
      signal.data[to.adjust.pos,i] <- max.value
      to.adjust.pos <- which(signal.data[,i] <= min.value)
      signal.data[to.adjust.pos,i] <- min.value
    }
  }
  #Remove those values selected
  if (length(to.remove.pos.all) > 0){
    signal.data <- signal.data[-to.remove.pos.all,]
  }
  result <- list(signal.data, to.remove.pos.all)
  return(result)
}

#Function to label the cells if those are cancer-cells or not according to the expression of CK7, ECadherin, Vimentin and aSMA
gaussian.mix.cancer.cell.detection <- function(dataset, cancer.channels = c("CK7", "ECadherin"), stromal.channels = c("Vimentin", "aSMA"), 
                                               patientIDstable=cores.patient){
  data.filtered.patient.scaled <- dataset
  data.filtered.patient.scaled$is.Cancer.cell <- NA
  #Scaling by patient the selected channels (channels norm)
  for (pt in unique(patientIDstable$patient)){
    #pt <- unique(patientIDstable$patient)[1]
    cores.pt <- as.numeric(gsub("core","",patientIDstable[patientIDstable$patient==pt, "cycif.core.id"]))
    
    #Scaling data signal from stromal channel markers (Vimentin and aSMA)
    stromal.channels.scaled <- data.filtered.patient.scaled %>%
      dplyr::filter(CoreId %in% cores.pt) %>%
      dplyr::select(stromal.channels)
    data.trim <- z.trimming(stromal.channels.scaled)
    stromal.channels.scaled <- data.trim[[1]]

    #Getting the median and sd values for Vimentin and aSMA per patient
    cut.value.vimentin <- mean(stromal.channels.scaled$Vimentin) + sd(stromal.channels.scaled$Vimentin)
    cut.value.aSMA <- mean(stromal.channels.scaled$aSMA) + sd(stromal.channels.scaled$aSMA)
    
    #Scaling data signal from cancer channel markers (CK7 and ECadherin)
    cancer.channels.scaled <- data.filtered.patient.scaled %>%
      dplyr::filter(CoreId %in% cores.pt) %>%
      dplyr::select(cancer.channels) %>%
      mutate_at(cancer.channels, ~(scale(.) %>% as.vector))
    data.trim <- z.trimming(cancer.channels.scaled)
    cancer.channels.scaled <- data.trim[[1]]
    cancer.channels.scaled2 <- apply(cancer.channels.scaled, 2, min_max_norm)
    
    #Fitting mixture of two Gaussians
    mod1 <- Mclust(cancer.channels.scaled2, G = 2, modelName = c("VVV","VVI"))
    summaryMclust <- summary(mod1, parameters = TRUE)
    cancer.lab <- which.max(summaryMclust$mean[1,]) #Cancer cells have the max of CK7 or ECadherin
    non.cancer.lab <- which.min(summaryMclust$mean[1,]) #Cancer cells have the max of CK7 or ECadherin
    channel.nanmes.Mclust <- row.names(summaryMclust$mean)
    
    #Cells which intensity of CK7 and ECadherin are below the mean intensity of non.cancer.cells, are labeled as non.cancer.cell
    adittional.non.cancer.rows <- which(cancer.channels.scaled2[,channel.nanmes.Mclust[1]] <= summaryMclust$mean[1,non.cancer.lab] &
                                                cancer.channels.scaled2[,channel.nanmes.Mclust[2]] <= summaryMclust$mean[2,non.cancer.lab])
    mod1$classification[adittional.non.cancer.rows] <- non.cancer.lab
    
    #Check if the signal intensity of Vimentin or aSMA is high on the classified cancer cells
    #If Vimentin or aSMA is above one sd above the median of all cells by patient and also lower than the mean cancer expression of CK7 and ECadherin in cancer cells
    #Then those cells would be classified as non-cancer
    false.positive.cancer.rows <-    which(mod1$classification == cancer.lab &
                                      (stromal.channels.scaled$Vimentin >= cut.value.vimentin | 
                                       stromal.channels.scaled$aSMA >= cut.value.aSMA) &
                                       cancer.channels.scaled2[,channel.nanmes.Mclust[1]] < summaryMclust$mean[1,cancer.lab] &
                                       cancer.channels.scaled2[,channel.nanmes.Mclust[2]] < summaryMclust$mean[2,cancer.lab])
    
    mod1$classification[false.positive.cancer.rows] <- non.cancer.lab
    
    data.filtered.patient.scaled[data.filtered.patient.scaled$CoreId %in% cores.pt, "is.Cancer.cell"] <- ifelse(mod1$classification == cancer.lab, "YES","NO")
    
    png(paste0(outout.folder.name,pt,"_gaussians-class.png"), width = 640, height = 560)
    plot(mod1, what = "classification", main=pt)
    text(0.8,0.7,pt)
    dev.off()
  }
  data.filtered.patient.scaled <- data.filtered.patient.scaled[,c("Row_number","is.Cancer.cell")]
  return(data.filtered.patient.scaled)
}


#Function to gate the cells if those are cancer-cells or not
gating.cancer.cell.detection <- function(dataset, cancer.channels = c("CK7", "ECadherin"), stromal.channels = c("Vimentin", "aSMA"), 
                                               patientIDstable=cores.patient, gating.file=to.gate.cancer.cells){
  data.filtered.patient.scaled <- dataset
  data.filtered.patient.scaled$is.Cancer.cell <- NA
  #Scaling by patient the selected channels (channels norm)
  for (pt in unique(patientIDstable$patient)){
    cores.pt <- as.numeric(gsub("core","",patientIDstable[patientIDstable$patient==pt, "cycif.core.id"]))
    gates.pt <- gating.file[gating.file$Patient %in% pt,]
    
    #Scaling data signal from stromal channel markers (Vimentin and aSMA)
    stromal.channels.scaled <- data.filtered.patient.scaled %>%
      dplyr::filter(CoreId %in% cores.pt) %>%
      dplyr::select(stromal.channels)
    data.trim <- z.trimming(stromal.channels.scaled)
    stromal.channels.scaled <- data.trim[[1]]
    
    #Getting the median and sd values for Vimentin and aSMA per patient
    cut.value.vimentin <- mean(stromal.channels.scaled$Vimentin) + sd(stromal.channels.scaled$Vimentin)
    cut.value.aSMA <- mean(stromal.channels.scaled$aSMA) + sd(stromal.channels.scaled$aSMA)
    
    #Scaling data signal from cancer channel markers (CK7 and ECadherin)
    cancer.channels.scaled <- data.filtered.patient.scaled %>%
      dplyr::filter(CoreId %in% cores.pt) %>%
      dplyr::select(cancer.channels) %>%
      mutate_at(cancer.channels, ~(scale(.) %>% as.vector))
    data.trim <- z.trimming(cancer.channels.scaled)
    cancer.channels.scaled <- data.trim[[1]]
    cancer.channels.scaled2 <- apply(cancer.channels.scaled, 2, min_max_norm)
    
    #Gating the cells according the the CK7 and ECadherin signal intensity
    cancer.cell.status <-  ifelse(cancer.channels.scaled2[,"CK7"] <= gates.pt$CK7 & cancer.channels.scaled2[,"ECadherin"] <= gates.pt$Ecadherin, "NO","YES")
    
    #Getting the median values for CK7 and ECadherin in the cancer cells
    CK7.med.cancer <- median(cancer.channels.scaled2[cancer.cell.status == "YES","CK7"])
    ECadherin.med.cancer <- median(cancer.channels.scaled2[cancer.cell.status == "YES","ECadherin"])
      
    #Check if the signal intensity of Vimentin or aSMA is high on the classified cancer cells
    #If Vimentin or aSMA is above one sd above the median of all cells by patient and also lower than the median cancer expression of CK7 and ECadherin in cancer cells
    #Then those cells would be classified as non-cancer
    false.positive.cancer.rows <- which(cancer.cell.status == "YES" &
                                             (stromal.channels.scaled$Vimentin >= cut.value.vimentin | 
                                                stromal.channels.scaled$aSMA >= cut.value.aSMA) &
                                             cancer.channels.scaled2[,"CK7"] < CK7.med.cancer &
                                             cancer.channels.scaled2[,"ECadherin"] < ECadherin.med.cancer)
    cancer.cell.status[false.positive.cancer.rows] <- "NO"
    
    data.filtered.patient.scaled[data.filtered.patient.scaled$CoreId %in% cores.pt, "is.Cancer.cell"] <- cancer.cell.status
    
    color.for.plor <- ifelse(cancer.cell.status == "YES", "red", "royalblue")
    pch.for.plor <- ifelse(cancer.cell.status == "YES", 0, 16)
    
    png(paste0(outout.folder.name,pt,"_gaussians-class.png"), width = 640, height = 560)
    plot(cancer.channels.scaled2[,"CK7"],  cancer.channels.scaled2[,"ECadherin"], main=pt, xlab="CK7", ylab="ECadherin",
         pch=pch.for.plor, col=color.for.plor)
    text(0.8,0.7,pt)
    dev.off()
  }
  data.filtered.patient.scaled <- data.filtered.patient.scaled[,c("Row_number","is.Cancer.cell")]
  return(data.filtered.patient.scaled)
}
