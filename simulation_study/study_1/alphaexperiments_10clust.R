library(tidyverse)
library(mclust)
library(VICatMix)

alpha <- c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 5) #values of alpha to test
set.seed(13445)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)
generation <- generateSampleDataBin(1000, 10, rep(0.1, 10), 100, 0)
data <- generation[[1]]
truelabels <- generation[[2]]
alphaexp1 <- foreach(i = 1:10) %dorng% {
  #Simulated data n = 1000, clusters = 10, variables = 20
  run <- lapply(X = alpha, FUN = runVICatMix, data = data, K = 30, tol = 0.00005)
  finalELBO <- data.frame(alpha = alpha, ELBO = NA, Clusters = NA)
  for (a in 1:length(alpha)){
    finalELBO$ELBO[a] <- run[[a]]$ELBO[length(run[[a]]$ELBO)]
    finalELBO$Clusters[a] <- run[[a]]$Cl[length(run[[a]]$Cl)]
    finalELBO$ARI[a] <- adjustedRandIndex(truelabels$Cluster, run[[a]]$model$labels)
  }
  finalELBO
}

save(alphaexp1, file="alphaexp1.RData")

