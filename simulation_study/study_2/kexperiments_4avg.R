#Changing initial value of K in VICatMix-Avg (Section S5.2.2) - 4 true clusters
library(tidyverse)
library(mclust)
library(VICatMix)

K_vals <- c(3, 4, 5, 6, 7, 8, 9, 10) #values of alpha to test
set.seed(711)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(25)
Kexp1_avg <- list()
for (j in 1:10){
  generation <- generateSampleDataBin(1000, 4, c(0.1, 0.2, 0.3, 0.4), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  Kexp1_avgindiv <- list()
  for(i in 1:5){
    #Simulated data n = 1000, clusters = 4, variables = 100
    run <- lapply(X = K_vals, FUN = runVICatMixAvg, data = data, alpha = 0.01, parallel = TRUE)
    finalELBO <- data.frame(K = K_vals, Clusters = NA)
    for (a in 1:length(K_vals)){
      finalELBO$Clusters[a] <- length(unique(run[[a]]$labels_avg))
      finalELBO$ARI[a] <- adjustedRandIndex(truelabels$Cluster, run[[a]]$labels_avg)
    }
    Kexp1_avgindiv[[i]] <- finalELBO
  }
  Kexp1_avg[[j]] <- Kexp1_avgindiv
}

save(Kexp1_avg, file="Kexp1_avg.RData")
