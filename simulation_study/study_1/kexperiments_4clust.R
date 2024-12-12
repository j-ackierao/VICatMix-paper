library(tidyverse)
library(mclust)
library(VICatMix)

K_vals <- c(3, 4, 5, 6, 7, 8, 9, 10) #values of alpha to test
set.seed(13445)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)
Kexp1 <- list()
for (j in 1:10){
  generation <- generateSampleDataBin(1000, 4, c(0.1, 0.2, 0.3, 0.4), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  Kexp1_indiv <- foreach(i = 1:10) %dorng% {
    #Simulated data n = 1000, clusters = 4, variables = 20
    run <- lapply(X = K_vals, FUN = runVICatMix, data = data, alpha = 0.01)
    finalELBO <- data.frame(K = K_vals, Clusters = NA)
    for (a in 1:length(K_vals)){
      finalELBO$Clusters[a] <- run[[a]]$Cl[length(run[[a]]$Cl)]
      finalELBO$ARI[a] <- adjustedRandIndex(truelabels$Cluster, run[[a]]$model$labels)
    }
    finalELBO
  }
  Kexp1[[j]] <- Kexp1_indiv
}

save(Kexp1, file="Kexp1.RData")

