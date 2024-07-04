#Comparing run times with increasing n 

library(tidyverse)
library(mclust)
source("Variational Mixture Model.R")
source("GenerateSampleData.R")

set.seed(355519)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)

n_list <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
n_runtimes <- foreach(i = 1:10) %dorng% {
  runtimes_n <- list()
  for (j in 1:9){
    generation <- GenerateSampleData(n_list[j], 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), 100, 0)
    data <- generation[[1]]
    truelabels <- generation[[2]]
    start.time <- Sys.time()
    clustering <- mixturemodel(data, 20, 0.05, 2000, 0.00005)
    end.time <- Sys.time()
    timetaken <- difftime(end.time, start.time, unit = "mins")
    modelresult <- data.frame(n_val = n_list[j], ELBO = NA, Clusters = NA, ARI = NA, Time = NA)
    modelresult$ELBO[1] <- clustering$ELBO[length(clustering$ELBO)]
    modelresult$Clusters[1] <- clustering$Cl[length(clustering$Cl)]
    modelresult$ARI[1] <- adjustedRandIndex(truelabels$Cluster, clustering$model$labels)
    modelresult$Time[1] <- timetaken
    runtimes_n[[j]] <- modelresult
  }
  runtimes_n
}
