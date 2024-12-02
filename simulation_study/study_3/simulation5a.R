#Simulation 2 (variable selection - 50/50, n = 2000)
#BHC is in a different file (simulation5b) due to memory usage

library(tidyverse)
library(mclust)
library(BHC)
library(PReMiuM)
library(factoextra)
library(stats)
library(mcclust)
library(mcclust.ext)
library(VICatMix)
set.seed(1931253)

#Need function to generate BHC clusters


library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)
sim2varsel <- foreach(i = 1:10) %dorng% {
  generation <- generateSampleDataBin(2000, 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), 50, 50)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  itemLabels <- rownames(data)
  colnames(data) <- paste0("Covariate",1:100)
  final_result <- list()
  results_list <- list()
  model_compare <- data.frame(model = c("VICatMix", "VICatMix-Avg", "VICatMixVarSel", "VICatMixVarSel-Avg", "PReMiuM-NoVarSel", "PReMiuM-BinVarSel", "PReMiuM-CtsVarSel"), 
                              ARI = NA, Clusters = NA, Time = NA, Rel = NA, Irrel = NA, Rel5 = NA, Irrel5 = NA, RelMean = NA, IrrelMean = NA, RelMean5 = NA, IrrelMean5 = NA)
  #Try every model
  #VICatMix
  vicatmix_labels <- list() #save labels 
  vicatmix_ELBO <- c()
  vicatmix_time <- c()
  vicatmix_clust <- c()
  for (j in 1:25){
    start.time <- Sys.time()
    vicatmix <- runVICatMix(data, 20, 0.05, tol = 0.000005)
    end.time <- Sys.time()
    vicatmix_time[j] <- difftime(end.time, start.time, unit = "secs")
    vicatmix_labels[[j]] <- vicatmix$model$labels 
    vicatmix_clust[j] <- length(unique(vicatmix$model$labels ))
    vicatmix_ELBO[j] <- vicatmix$ELBO[length(vicatmix$ELBO)]
  }
  p1 <- t(matrix(unlist(vicatmix_labels), 2000, 25))
  psm <- comp.psm(p1)
  VIcomp <- minVI(psm, method = 'comp',max.k = 20)$cl
  bestVI <- which.max(vicatmix_ELBO) 
  results_list$VICatMixLabels <- vicatmix_labels[[bestVI]]
  model_compare$ARI[1] <- adjustedRandIndex(vicatmix_labels[[bestVI]], truelabels$Cluster)
  model_compare$Time[1] <- vicatmix_time[[bestVI]]
  model_compare$Clusters[1] <- length(unique(vicatmix_labels[[bestVI]]))
  results_list$VICatMixAvgLabels <- VIcomp
  model_compare$ARI[2] <- adjustedRandIndex(VIcomp, truelabels$Cluster)
  model_compare$Time[2] <- NA
  model_compare$Clusters[2] <- length(unique(VIcomp))
  
  vicatmix_labels$ELBO <- vicatmix_ELBO
  vicatmix_labels$Clust <- vicatmix_clust
  vicatmix_labels$time <- vicatmix_time
  
  #VICatMixVarSel
  vicatmixvs_labels <- list() #save labels 
  vicatmixvs_ELBO <- c()
  vicatmixvs_time <- c()
  vicatmixvs_clust <- c()
  vicatmixvs_vars <- list() #save selected variables
  for (j in 1:25){
    start.time <- Sys.time()
    vicatmixvs <- runVICatMixVarSel(data, 20, 0.05, tol = 0.000005)
    end.time <- Sys.time()
    vicatmixvs_time[j] <- difftime(end.time, start.time, unit = "secs")
    vicatmixvs_labels[[j]] <- vicatmixvs$model$labels 
    vicatmixvs_clust[j] <- length(unique(vicatmixvs$model$labels ))
    vicatmixvs_ELBO[j] <- vicatmixvs$ELBO[length(vicatmixvs$ELBO)]
    vicatmixvs_vars[[j]] <- vicatmixvs$model$c
  }
  p1 <- t(matrix(unlist(vicatmixvs_labels), 2000, 25))
  psm <- comp.psm(p1)
  VIcompvs <- minVI(psm, method = 'comp',max.k = 20)$cl
  bestVI <- which.max(vicatmixvs_ELBO) 
  results_list$VICatMixVarSelLabels <- vicatmixvs_labels[[bestVI]]
  model_compare$ARI[3] <- adjustedRandIndex(vicatmixvs_labels[[bestVI]], truelabels$Cluster)
  model_compare$Time[3] <- vicatmixvs_time[[bestVI]]
  model_compare$Clusters[3] <- length(unique(vicatmixvs_labels[[bestVI]]))
  model_compare$Rel[3] <- sum(vicatmixvs_vars[[bestVI]][1:50] >= 0.95)
  model_compare$Irrel[3] <- sum(vicatmixvs_vars[[bestVI]][51:100] < 0.95)
  results_list$VICatMixVarSelAvgLabels <- VIcompvs
  model_compare$ARI[4] <- adjustedRandIndex(VIcompvs, truelabels$Cluster)
  model_compare$Time[4] <- NA
  model_compare$Clusters[4] <- length(unique(VIcompvs))
  vivarss <- t(matrix(unlist(vicatmixvs_vars), 100, 25))
  vivarss[vivarss > 0.5] <- 1
  vivarss[vivarss <= 0.5] <- 0
  result <- vector(mode = 'numeric', length = 100)
  for (j in 1:100){
    result[j] <- sum(vivarss[,j]) / 25
  }
  model_compare$Rel[4] <- sum(result[1:50] >= 0.95)
  model_compare$Irrel[4] <- sum(result[51:100] < 0.95)
  
  vicatmixvs_labels$ELBO <- vicatmixvs_ELBO
  vicatmixvs_labels$Clust <- vicatmixvs_clust
  vicatmixvs_labels$time <- vicatmixvs_time
  vicatmixvs_labels$vars <- vicatmixvs_vars
  
  #PReMiuM
  
  covNames  <- colnames(data)
  outcome   <- c(seq(0,0,length = ceiling(nrow(data)/2)), seq(1,1,length = floor(nrow(data)/2)))
  outdata <- cbind(outcome, data) # This is the data file (includes response and covariates)
  
  # Put data and modelling options into an "inputs" data frame
  inputs                  <- generateSampleDataFile(clusSummaryBernoulliDiscrete()) # This is just to initialise the "inputs" variable, so that it has all of the right fieldnames
  inputs$inputData        <- as.data.frame(outdata)   #Put our data in the "inputData" slot
  inputs$covNames         <- covNames
  inputs$nCovariates      <- length(covNames)
  inputs$fixedEffectNames <- NULL
  
  start.time <- Sys.time()
  premiumModel <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("sim_output/varselsim2_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE)
  dissimObj <-calcDissimilarityMatrix(premiumModel)
  clusObj   <-calcOptimalClustering(dissimObj, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMlabels <- clusObj$clustering
  model_compare$ARI[5] <- adjustedRandIndex(truelabels$Cluster, clusObj$clustering)
  model_compare$Time[5] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[5] <- length(unique(clusObj$clustering))
  
  #PReMiuM - BinVarSel
  
  start.time <- Sys.time()
  premiumModel2 <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("output/varselsim2bin_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE, varSelectType = 'BinaryCluster', seed=sample.int(2^32-1, 1))
  dissimObj2 <-calcDissimilarityMatrix(premiumModel2)
  clusObj2   <-calcOptimalClustering(dissimObj2, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMBinlabels <- clusObj2$clustering
  model_compare$ARI[6] <- adjustedRandIndex(truelabels$Cluster, clusObj2$clustering)
  model_compare$Time[6] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[6] <- length(unique(clusObj2$clustering))
  BinRhoMedian <- summariseVarSelectRho(premiumModel2)$rhoMedian
  BinRhoMean <- summariseVarSelectRho(premiumModel2)$rhoMean
  model_compare$Rel[6] <- sum(BinRhoMedian[1:50] >= 0.95)
  model_compare$Irrel[6] <- sum(BinRhoMedian[51:100] < 0.95)
  model_compare$Rel5[6] <- sum(BinRhoMedian[1:50] >= 0.5)
  model_compare$Irrel5[6] <- sum(BinRhoMedian[51:100] < 0.5)
  model_compare$RelMean[6] <- sum(BinRhoMean[1:50] >= 0.95)
  model_compare$IrrelMean[6] <- sum(BinRhoMean[51:100] < 0.95)
  model_compare$RelMean5[6] <- sum(BinRhoMean[1:50] >= 0.5)
  model_compare$IrrelMean5[6] <- sum(BinRhoMean[51:100] < 0.5)
  
  #PReMiuM - CtsVarSel
  
  start.time <- Sys.time()
  premiumModel3 <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("output/varselsim2cts_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE, varSelectType = 'Continuous', seed=sample.int(2^32-1, 1))
  dissimObj3 <-calcDissimilarityMatrix(premiumModel3)
  clusObj3   <-calcOptimalClustering(dissimObj3, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMCtslabels <- clusObj2$clustering
  model_compare$ARI[7] <- adjustedRandIndex(truelabels$Cluster, clusObj3$clustering)
  model_compare$Time[7] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[7] <- length(unique(clusObj3$clustering))
  CtsRhoMedian <- summariseVarSelectRho(premiumModel3)$rhoMedian
  CtsRhoMean <- summariseVarSelectRho(premiumModel3)$rhoMean
  model_compare$Rel[7] <- sum(CtsRhoMedian[1:50] >= 0.95)
  model_compare$Irrel[7] <- sum(CtsRhoMedian[51:100] < 0.95)
  model_compare$Rel5[7] <- sum(CtsRhoMedian[1:50] >= 0.5)
  model_compare$Irrel5[7] <- sum(CtsRhoMedian[51:100] < 0.5)
  model_compare$RelMean[7] <- sum(CtsRhoMean[1:50] >= 0.95)
  model_compare$IrrelMean[7] <- sum(CtsRhoMean[51:100] < 0.95)
  model_compare$RelMean5[7] <- sum(CtsRhoMean[1:50] >= 0.5)
  model_compare$IrrelMean5[7] <- sum(CtsRhoMean[51:100] < 0.5)
  
  
  final_result$model_compare <- model_compare
  final_result$model_labels <- results_list
  final_result$VIres <- vicatmix_labels
  final_result$VIvarselres <- vicatmixvs_labels
  final_result$generation <- generation
  final_result
}

save(sim2varsel, file="/home/jr951/VariationalMixtures/sim2varsel.RData")