library(tidyverse)
library(mclust)
library(VICatMix)

K_vals <- c(5, 10, 15, 20, 25, 30) #values of alpha to test
set.seed(13445)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)
Kexp2 <- list()
for (j in 1:10){
  generation <- generateSampleDataBin(1000, 10, rep(0.1, 10), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  Kexp2_indiv <- foreach(i = 1:10) %dorng% {
    #Simulated data n = 1000, clusters = 4, variables = 20
    run <- lapply(X = K_vals, FUN = runVICatMix, data = data, alpha = 0.01)
    finalELBO <- data.frame(K = K_vals, Clusters = NA)
    for (a in 1:length(K_vals)){
      finalELBO$Clusters[a] <- run[[a]]$Cl[length(run[[a]]$Cl)]
      finalELBO$ARI[a] <- adjustedRandIndex(truelabels$Cluster, run[[a]]$model$labels)
    }
    finalELBO
  }
  Kexp2[[j]] <- Kexp2_indiv
}

save(Kexp2, file="Kexp2.RData")

