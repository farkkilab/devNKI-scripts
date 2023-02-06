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



# dataset <- data.filtered
# patientIDstable = cores.patient
# other.chanels.norm = chanels.global.gates.non.cancer
# chanels.norm = c("CK7", "ECadherin")

#Normalization of CK7 and ECadherin by patient
normalize.by.patient <- function(dataset, chanels.norm = c("CK7", "ECadherin"), other.chanels.norm=NULL,
                                   patientIDstable=cores.patient){
  data.filtered.patient.scaled <- dataset
  data.filtered.patient.scaled$is.Cancer.cell <- NA
  #Scaling by patient the selected channels (channels norm)
  for (pt in unique(patientIDstable$patient)){
    cores.pt <- as.numeric(gsub("core","",patientIDstable[patientIDstable$patient==pt, "cycif.core.id"]))
    selected.channels.scaled <- data.filtered.patient.scaled %>%
                    dplyr::filter(CoreId %in% cores.pt) %>%
                    dplyr::select(chanels.norm) %>%
                    mutate_at(chanels.norm, ~(scale(.) %>% as.vector))
    data.trim <- z.trimming(selected.channels.scaled)
    selected.channels.scaled <- data.trim[[1]]
    selected.channels.scaled2 <- BBmisc::normalize(selected.channels.scaled, method='range')
    
    #Fitting mixture of two gaussians
    mod1 <- Mclust(selected.channels.scaled2, G = 2, modelName = c("VVV","VVI"))
    summaryMclust <- summary(mod1, parameters = TRUE)
    cancer.lab <- which.max(summaryMclust$mean[1,]) #Cancer cells have the max of CK7 or ECadherin
    non.cancer.lab <- which.min(summaryMclust$mean[1,]) #Cancer cells have the max of CK7 or ECadherin
    channel.nanmes.Mclust <- row.names(summaryMclust$mean)
    
    #Cells which intensity of CK7 and ECadherin are below the mean intensity of non.cancer.cells, are labeled as non.cancer.cell
    adittional.non.cancer.rows <- which(selected.channels.scaled2[,channel.nanmes.Mclust[1]] <= summaryMclust$mean[1,non.cancer.lab] & selected.channels.scaled2[,channel.nanmes.Mclust[2]] <= summaryMclust$mean[2,non.cancer.lab])
    mod1$classification[adittional.non.cancer.rows] <- non.cancer.lab
    
    data.filtered.patient.scaled[data.filtered.patient.scaled$CoreId %in% cores.pt, "is.Cancer.cell"] <- ifelse(mod1$classification == cancer.lab, "YES","NO")
    png(paste0(outout.folder.name,pt,"_gaussians-class.png"), width = 640, height = 560)
    plot(mod1, what = "classification", main=pt)
    text(0.8,0.7,pt)
    dev.off()
  }
  #Scaling the rest of the channels
  # data.filtered.patient.scaled[,other.chanels.norm] <- scale(data.filtered.patient.scaled[,other.chanels.norm])
  # data.filtered.patient.scaled[,other.chanels.norm] <- BBmisc::normalize(data.filtered.patient.scaled[,other.chanels.norm], method='range')
  return(data.filtered.patient.scaled)
}
