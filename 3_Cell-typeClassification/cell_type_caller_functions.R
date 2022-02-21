###A function to truncate the dataset, using the 99.9 percentile
#Other commented lines in the script were added to make the truncation using percentiles and ignores cells with a Z score so high
z.trimming <- function(signal.data, z.cutoff = 10, z.for.max.intensity = 2){
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
      max.value <- quantile(signal.data[,i], 0.999)
      min.value <- quantile(signal.data[,i], 0.001)
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