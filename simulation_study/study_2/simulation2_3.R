#Summarisation of multiple runs of variational model: Tests with normal PS<

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("Variational Mixture Model.R")
source("GenerateSampleData.R")

set.seed(992200)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)

psmtests2000 <- foreach(i = 1:10) %dorng% {
  generation <- GenerateSampleData(2000, 20, rep(c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), 2), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  resultforpsm <- list()
  for (j in 1:30){
    mix <- mixturemodel(data, 40, 0.05, 2000, 0.00005)
    resultforpsm[[j]] <- mix$model$labels
  }
  run_nums <- c(5, 10, 15, 20, 25, 30)
  all_result <- list()
  for (k in 1:6){
    resultforpsm2 <- resultforpsm[1:run_nums[k]]
    p1 <- t(matrix(unlist(resultforpsm2), 2000, run_nums[k]))
    psm <- comp.psm(p1)
    start.time <- Sys.time()
    medv <- medv(psm)
    end.time <- Sys.time()
    medv.time <- difftime(end.time, start.time, unit = "mins")
    start.time <- Sys.time()
    VIcomp <- minVI(psm, method = 'comp',max.k = 40)$cl
    end.time <- Sys.time()
    VIcomp.time <- difftime(end.time, start.time, unit = "mins")
    start.time <- Sys.time()
    VIavg <- minVI(psm, method = 'avg',max.k = 40)$cl
    end.time <- Sys.time()
    VIavg.time <- difftime(end.time, start.time, unit = "mins")
    summ_result <- data.frame(model = c("medv", "VIcomp", "VIavg"), run_num = run_nums[k], Clusters = NA, ARI = NA, Time = NA)
    summ_result$ARI[1] <- adjustedRandIndex(medv, truelabels$Cluster)
    summ_result$ARI[2] <- adjustedRandIndex(VIcomp, truelabels$Cluster)
    summ_result$ARI[3] <- adjustedRandIndex(VIavg, truelabels$Cluster)
    summ_result$Clusters[1] <- length(unique(medv))
    summ_result$Clusters[2] <- length(unique(VIcomp))
    summ_result$Clusters[3] <- length(unique(VIavg))
    summ_result$Time[1] <- medv.time
    summ_result$Time[2] <- VIavg.time
    summ_result$Time[3] <- VIcomp.time
    all_result[[k]] <- summ_result
  }
  
  og_result <- data.frame(run_num = c(1:30), Clusters = NA, ARI = NA)
  for (j in 1:30){
    og_result$Clusters[j] <- length(unique(resultforpsm[[j]]))
    og_result$ARI[j] <- adjustedRandIndex(resultforpsm[[j]], truelabels$Cluster)
  }
  all_result[[7]] <- og_result
  all_result
}

save(psmtests2000, file="/home/jr951/VariationalMixtures/psmtests2000.RData")


