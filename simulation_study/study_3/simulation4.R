#Simulation 1 (variable selection - 75/25)

library(tidyverse)
library(mclust)
library(BHC)
library(PReMiuM)
library(flexmix)
library(factoextra)
library(stats)
library(mcclust)
library(mcclust.ext)
library(VICatMix)
set.seed(2321218)

#Need function to generate BHC clusters
getBHCclusters <- function (dendro, rowLabels) 
{
  ##for ease, we use discrete height labels here
  ##this hardwires the logEvidence threshold at zero, which is the point
  ##where merge and not-merge hypotheses are equally likely  
  WhereToCut <- function(n) {
    attr(n, "height") <- 1
    if (!is.leaf(n)) {
      attr(n, "height") <- 2
      if (attr(n, "logEvidence") < 0) 
        attr(n, "height") <- 3
    }
    n
  }
  dendro <- dendrapply(dendro, WhereToCut) #Apply function to all nodes of a dendogram
  #ie process dendrogram nodes recursively
  cutDendro <- cut(dendro, 2) #cuts dendrogram at height 2
  #cut the dendrogram
  nClusters <- length(cutDendro$lower) #number of clusters ie. number of leaves after cutting above the leaves
  nTotalLabels <- length(labels(dendro)) #number of leaves ie. number of data points
  allLabels    <- vector(mode = "character", length = nTotalLabels)
  allClusters <- vector(mode = "integer", length = nTotalLabels)
  
  myCounter <- 1
  
  for (i in 1:nClusters) {
    #cut current dendrogram
    currentCluster <- cutDendro$lower[[i]] #takes ith cluster
    currentLabels <- labels(currentCluster) #labels in the cluster
    nLabels <- length(currentLabels) #number of items in the clsuter
    
    for (j in 1:nLabels) {
      allLabels[myCounter]   <- currentLabels[j] #label of one of the datapoints
      allClusters[myCounter] <- i #which cluster we are in
      myCounter <- myCounter + 1 #iterate over every datapoint one by one so we insert values into allLabels one by one
    }
  }
  
  clusters        <- allClusters
  names(clusters) <- allLabels
  clusters <- clusters[rowLabels]
  clustersdf <- as.data.frame(clusters)
  
  return(clustersdf)
}

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(10)
sim1varsel <- foreach(i = 1:10) %dorng% {
  generation <- generateSampleDataBin(1000, 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), 75, 25)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  itemLabels <- rownames(data)
  colnames(data) <- paste0("Covariate",1:100)
  final_result <- list()
  results_list <- list()
  model_compare <- data.frame(model = c("VICatMix", "VICatMix-Avg", "VICatMixVarSel", "VICatMixVarSel-Avg", "PReMiuM-NoVarSel", "PReMiuM-BinVarSel", "PReMiuM-CtsVarSel", "FlexMixICL", "BHC"), 
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
  p1 <- t(matrix(unlist(vicatmix_labels), 1000, 25))
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
  p1 <- t(matrix(unlist(vicatmixvs_labels), 1000, 25))
  psm <- comp.psm(p1)
  VIcompvs <- minVI(psm, method = 'comp',max.k = 20)$cl
  bestVI <- which.max(vicatmixvs_ELBO) 
  results_list$VICatMixVarSelLabels <- vicatmixvs_labels[[bestVI]]
  model_compare$ARI[3] <- adjustedRandIndex(vicatmixvs_labels[[bestVI]], truelabels$Cluster)
  model_compare$Time[3] <- vicatmixvs_time[[bestVI]]
  model_compare$Clusters[3] <- length(unique(vicatmixvs_labels[[bestVI]]))
  model_compare$Rel[3] <- sum(vicatmixvs_vars[[bestVI]][1:75] >= 0.95)
  model_compare$Irrel[3] <- sum(vicatmixvs_vars[[bestVI]][76:100] < 0.95)
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
  model_compare$Rel[4] <- sum(result[1:75] >= 0.95)
  model_compare$Irrel[4] <- sum(result[76:100] < 0.95)
  
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
  premiumModel <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("sim_output/varselsim1_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE)
  dissimObj <-calcDissimilarityMatrix(premiumModel)
  clusObj   <-calcOptimalClustering(dissimObj, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMlabels <- clusObj$clustering
  model_compare$ARI[5] <- adjustedRandIndex(truelabels$Cluster, clusObj$clustering)
  model_compare$Time[5] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[5] <- length(unique(clusObj$clustering))
  
  #PReMiuM - BinVarSel
  
  start.time <- Sys.time()
  premiumModel2 <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("output/varselsim1bin_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE, varSelectType = 'BinaryCluster', seed=sample.int(2^32-1, 1))
  dissimObj2 <-calcDissimilarityMatrix(premiumModel2)
  clusObj2   <-calcOptimalClustering(dissimObj2, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMBinlabels <- clusObj2$clustering
  model_compare$ARI[6] <- adjustedRandIndex(truelabels$Cluster, clusObj2$clustering)
  model_compare$Time[6] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[6] <- length(unique(clusObj2$clustering))
  BinRhoMedian <- summariseVarSelectRho(premiumModel2)$rhoMedian
  BinRhoMean <- summariseVarSelectRho(premiumModel2)$rhoMean
  model_compare$Rel[6] <- sum(BinRhoMedian[1:75] >= 0.95)
  model_compare$Irrel[6] <- sum(BinRhoMedian[76:100] < 0.95)
  model_compare$Rel5[6] <- sum(BinRhoMedian[1:75] >= 0.5)
  model_compare$Irrel5[6] <- sum(BinRhoMedian[76:100] < 0.5)
  model_compare$RelMean[6] <- sum(BinRhoMean[1:75] >= 0.95)
  model_compare$IrrelMean[6] <- sum(BinRhoMean[76:100] < 0.95)
  model_compare$RelMean5[6] <- sum(BinRhoMean[1:75] >= 0.5)
  model_compare$IrrelMean5[6] <- sum(BinRhoMean[76:100] < 0.5)
  
  #PReMiuM - CtsVarSel
  
  start.time <- Sys.time()
  premiumModel3 <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("output/varselsim1cts_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE, varSelectType = 'Continuous', seed=sample.int(2^32-1, 1))
  dissimObj3 <-calcDissimilarityMatrix(premiumModel3)
  clusObj3   <-calcOptimalClustering(dissimObj3, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMCtslabels <- clusObj3$clustering
  model_compare$ARI[7] <- adjustedRandIndex(truelabels$Cluster, clusObj3$clustering)
  model_compare$Time[7] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[7] <- length(unique(clusObj3$clustering))
  CtsRhoMedian <- summariseVarSelectRho(premiumModel3)$rhoMedian
  CtsRhoMean<- summariseVarSelectRho(premiumModel3)$rhoMean
  model_compare$Rel[7] <- sum(CtsRhoMedian[1:75] >= 0.95)
  model_compare$Irrel[7] <- sum(CtsRhoMedian[76:100] < 0.95)
  model_compare$Rel5[7] <- sum(CtsRhoMedian[1:75] >= 0.5)
  model_compare$Irrel5[7] <- sum(CtsRhoMedian[76:100] < 0.5)
  model_compare$RelMean[7] <- sum(CtsRhoMean[1:75] >= 0.95)
  model_compare$IrrelMean[7] <- sum(CtsRhoMean[76:100] < 0.95)
  model_compare$RelMean5[7] <- sum(CtsRhoMean[1:75] >= 0.5)
  model_compare$IrrelMean5[7] <- sum(CtsRhoMean[76:100] < 0.5)
  
  #BHC
  start.time <- Sys.time()
  bhcmodel <- bhc(data, itemLabels, verbose=TRUE)
  end.time <- Sys.time()
  bhcclust <- getBHCclusters(bhcmodel, itemLabels)$clusters
  results_list$BHCLabels <- bhcclust
  model_compare$ARI[9] <- adjustedRandIndex(bhcclust, truelabels$Cluster)
  model_compare$Clusters[9] <- length(unique(bhcclust))
  model_compare$Time[9] <- difftime(end.time, start.time, unit = "secs")
  
  #flexmix
  start.time <- Sys.time()
  flex <- initFlexmix(data ~ 1, k = 1:20, model = FLXMCmvbinary(), control = list(minprior = 0), nrep = 10)
  flexicl <- getModel(flex, which = "ICL")
  end.time <- Sys.time()
  flexiclclust <- clusters(flexicl)
  results_list$flexmixICL <- flexiclclust
  model_compare$ARI[8] <- adjustedRandIndex(flexiclclust, truelabels$Cluster)
  model_compare$Time[8] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[8] <- length(unique(flexiclclust))
  
  
  final_result$model_compare <- model_compare
  final_result$model_labels <- results_list
  final_result$VIres <- vicatmix_labels
  final_result$VIvarselres <- vicatmixvs_labels
  final_result
}

save(sim1varsel, file="/home/jr951/VariationalMixtures/sim1varsel.RData")