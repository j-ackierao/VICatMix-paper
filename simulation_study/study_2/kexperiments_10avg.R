#Changing initial value of K in VICatMix-Avg (Section S5.2.2) - 10 true clusters
library(tidyverse)
library(mclust)
library(VICatMix)

K_vals <- c(5, 10, 15, 20, 25, 30) #values of alpha to test
set.seed(711)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(25)
Kexp2_avg <- list()
for (j in 1:10){
  generation <- generateSampleDataBin(1000, 10, rep(0.1, 10), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  Kexp2_avgindiv <- list()
  for(i in 1:5){
    #Simulated data n = 1000, clusters = 10, variables = 100
    run <- lapply(X = K_vals, FUN = runVICatMixAvg, data = data, alpha = 0.01, parallel = TRUE)
    finalELBO <- data.frame(K = K_vals, Clusters = NA)
    for (a in 1:length(K_vals)){
      finalELBO$Clusters[a] <- length(unique(run[[a]]$labels_avg))
      finalELBO$ARI[a] <- adjustedRandIndex(truelabels$Cluster, run[[a]]$labels_avg)
    }
    Kexp2_avgindiv[[i]] <- finalELBO
  }
  Kexp2_avg[[j]] <- Kexp2_avgindiv
}

save(Kexp2_avg, file="Kexp2_avg.RData")
