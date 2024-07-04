#Look at correlation between ELBO and ARI

library(tidyverse)
library(mclust)
source("GenerateSampleData.R")
source("Variational Mixture Model.R")
set.seed(494)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(5)
ELBOARIcorr <- foreach(i = 1:5) %dorng% {
  generation <- GenerateSampleData(1000, 10, rep(0.1, 10), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  results_df <- data.frame(run = c(1:20), ELBO = NA, ARI = NA)
  for (j in 1:20){
    novarsel <- mixturemodel(data, 20, 0.05, 2000, 0.000005)
    results_df$ELBO[j] <- novarsel$ELBO[length(novarsel$ELBO)]
    results_df$ARI[j] <- adjustedRandIndex(novarsel$model$labels, truelabels$Cluster)
  }
  results_df
}

save(ELBOARIcorr, file="/home/jr951/VariationalMixtures/ELBOARIcorr.RData")

