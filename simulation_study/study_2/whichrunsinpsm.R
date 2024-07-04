#Testing using different runs - for supplement

library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)
source("Variational Mixture Model.R")
source("GenerateSampleData.R")

set.seed(13224)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)

whichrunsinpsm <- foreach(i = 1:10) %dorng% {
  #10 independently generated datasets
  generation <- GenerateSampleData(1000, 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  resultforpsm <- list()
  for (j in 1:50){
    #50 possible runs for each dataset
    mix <- mixturemodel(data, 30, 0.05, 2000, 0.00005)
    resultforpsm[[j]] <- mix$model$labels
  }
  all_result <- list()
  for (k in 1:10){
    #10 different combinations of 25 runs
    sample <- sample(1:50,25, replace=F) 
    resultforpsm2 <- resultforpsm[sample]
    p1 <- t(matrix(unlist(resultforpsm2), 1000, 25))
    psm <- comp.psm(p1)
    VIcomp <- minVI(psm, method = 'comp',max.k = 30)$cl
    summ_result <- data.frame(model = "VIcomp", Clusters = NA, ARI = NA, ogmean = NA, ogvar = NA, ogclustmean = NA, ogclustvar = NA)
    summ_result$ARI[1] <- adjustedRandIndex(VIcomp, truelabels$Cluster)
    summ_result$Clusters[1] <- length(unique(VIcomp))
    sampleARI <- c()
    sampleclust <- c()
    for (j in 1:25){
      #get sample mean and sample variance of each set of 25 runs
      whichrun <- sample[j]
      sampleclust[j] <- length(unique(resultforpsm[[whichrun]]))
      sampleARI[j] <- adjustedRandIndex(resultforpsm[[whichrun]], truelabels$Cluster)
    }
    summ_result$ogmean[1] <- mean(sampleARI)
    summ_result$ogvar[1] <- var(sampleARI)
    summ_result$ogclustmean[1] <- mean(sampleclust)
    summ_result$ogclustvar[1] <- var(sampleclust)
    all_result[[k]] <- summ_result
  }
  
  #Cluster numbers and ARIs for all 50 runs
  og_result <- data.frame(run_num = c(1:50), Clusters = NA, ARI = NA)
  for (j in 1:30){
    og_result$Clusters[j] <- length(unique(resultforpsm[[j]]))
    og_result$ARI[j] <- adjustedRandIndex(resultforpsm[[j]], truelabels$Cluster)
  }
  all_result[[11]] <- og_result
  all_result
}

save(whichrunsinpsm, file="/home/jr951/VariationalMixtures/whichrunsinpsm.RData")
