#Changing the value of a
library(tidyverse)
library(mclust)
source("VariationalMixtureModelVarSel.R")
source("GenerateSampleData.R")

cvals <- c(0.1, 0.5, 1, 2, 5) #values of a to test
set.seed(1621236)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)
cvals_sim <- foreach(i = 1:10) %dorng% {
  #Generate 10 independent datasets
  generation <- GenerateSampleData(1000, 10, rep(0.1, 10), 75, 25)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  results_list <- list()
  for(j in 1:10){
    #10 different initialisations on each dataset
    run <- lapply(X = cvals, FUN = mixturemodelvarsel, data = data, K = 20, alpha = 0.1, maxiter = 2000, tol = 0.00005)
    finalELBO <- data.frame(c = cvals, ELBO = NA, Clusters = NA, ARI = NA)
    for (a in 1:length(cvals)){
      finalELBO$ELBO[a] <- run[[a]]$ELBO[length(run[[a]]$ELBO)]
      finalELBO$Clusters[a] <- run[[a]]$Cl[length(run[[a]]$Cl)]
      finalELBO$ARI[a] <- adjustedRandIndex(truelabels$Cluster, run[[a]]$model$labels)
    }
    results_list[[j]] <- finalELBO
  }
  results_list
}

save(cvals_sim, file="/home/jr951/VariationalMixtures/cvals_sim.RData")