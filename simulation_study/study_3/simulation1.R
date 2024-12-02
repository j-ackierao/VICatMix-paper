#Simulation 1 (no variable selection)

library(tidyverse)
library(mclust)
library(BHC)
library(PReMiuM)
library(BayesBinMix)
library(flexmix)
library(factoextra)
library(stats)
library(mcclust)
library(mcclust.ext)
library(VICatMix)
set.seed(2321214)

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
sim1novarsel <- foreach(i = 1:10) %dorng% {
  
  generation <- generateSampleDataBin(1000, 10, rep(0.1, 10), 100, 0)
  data <- generation[[1]]
  truelabels <- generation[[2]]
  itemLabels <- rownames(data)
  colnames(data) <- paste0("Covariate",1:100)
  final_result <- list()
  results_list <- list()
  model_compare <- data.frame(model = c("VICatMix", "PReMiuM", "BHC", "BayesBinMix", "FlexMixICL", "FlexMixBIC", "hclust", "VICatMix-Avg"), 
                              ARI = NA, Clusters = NA, Time = NA)
  #Try every model
  #VICatMix
  #Manually run the model averaging so we can extract model with the best ELBO 
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
  model_compare$ARI[8] <- adjustedRandIndex(VIcomp, truelabels$Cluster)
  model_compare$Time[8] <- NA
  model_compare$Clusters[8] <- length(unique(VIcomp))
  
  vicatmix_labels$ELBO <- vicatmix_ELBO
  vicatmix_labels$Clust <- vicatmix_clust
  vicatmix_labels$time <- vicatmix_time
  
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
  premiumModel <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, data=inputs$inputData, nSweeps=2500, nClusInit=20, nBurn=1000, output=paste0("sim_output/novarselsim_", i), covNames = inputs$covNames, reportBurnIn = TRUE, excludeY = TRUE)
  dissimObj <-calcDissimilarityMatrix(premiumModel)
  clusObj   <-calcOptimalClustering(dissimObj, maxNClusters=20)
  end.time <- Sys.time()
  results_list$PReMlabels <- clusObj$clustering
  model_compare$ARI[2] <- adjustedRandIndex(truelabels$Cluster, clusObj$clustering)
  model_compare$Time[2] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[2] <- length(unique(clusObj$clustering))
  
  #BHC
  start.time <- Sys.time()
  bhcmodel <- bhc(data, itemLabels, verbose=TRUE)
  end.time <- Sys.time()
  bhcclust <- getBHCclusters(bhcmodel, itemLabels)$clusters
  results_list$BHCLabels <- bhcclust
  model_compare$ARI[3] <- adjustedRandIndex(bhcclust, truelabels$Cluster)
  model_compare$Clusters[3] <- length(unique(bhcclust))
  model_compare$Time[3] <- difftime(end.time, start.time, unit = "secs")
  
  #BayesBinMix
  start.time <- Sys.time()
  bayesbin <- coupledMetropolis(20, 6, seq(1, 0.6, length = 6), data, outPrefix = paste0("BayesBinSim_", i), ClusterPrior = "poisson", m = 2500, burn = 500)
  end.time <- Sys.time()
  bbclust <- bayesbin$clusterMembershipPerMethod$ECR
  results_list$BayesBinLabels <- bbclust
  model_compare$ARI[4] <- adjustedRandIndex(bbclust, truelabels$Cluster)
  model_compare$Time[4] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[4] <- length(unique(bbclust))
  
  bayesbin_labels <- bayesbin$clusterMembershipPerMethod #save labels 
  if(adjustedRandIndex(bayesbin_labels$STEPHENS, bayesbin_labels$ECR) != 1 |
     adjustedRandIndex(bayesbin_labels$STEPHENS, bayesbin_labels$ECR.ITERATIVE.1) != 1 |
     adjustedRandIndex(bayesbin_labels$ECR.ITERATIVE.1, bayesbin_labels$ECR) != 1){
    bayesbin_labels$DIFF <- "diff"
  }
  
  #flexmix
  start.time <- Sys.time()
  flex <- initFlexmix(data ~ 1, k = 1:20, model = FLXMCmvbinary(), control = list(minprior = 0), nrep = 10)
  flexicl <- getModel(flex, which = "ICL")
  end.time <- Sys.time()
  flexbic <- getModel(flex, which = "BIC")
  flexiclclust <- clusters(flexicl)
  flexbicclust <- clusters(flexbic)
  results_list$flexmixICL <- flexiclclust
  results_list$flexmixBIC <- flexbicclust
  model_compare$ARI[5] <- adjustedRandIndex(flexiclclust, truelabels$Cluster)
  model_compare$Time[5] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[5] <- length(unique(flexiclclust))
  model_compare$ARI[6] <- adjustedRandIndex(flexbicclust, truelabels$Cluster)
  model_compare$Time[6] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[6] <- length(unique(flexbicclust))
  
  #hclust
  start.time <- Sys.time()
  d <- dist(data)
  hier <- hclust(d)
  findsilhouette <- fviz_nbclust(data, FUN = hcut, method = 'silhouette')
  numclust <- which.max(findsilhouette$data$y)
  hierclust <- cutree(hier, k = numclust)
  end.time <- Sys.time()
  results_list$hclust <- hierclust
  model_compare$ARI[7] <- adjustedRandIndex(hierclust, truelabels$Cluster)
  model_compare$Time[7] <- difftime(end.time, start.time, unit = "secs")
  model_compare$Clusters[7] <- length(unique(hierclust))
  
  #hierarchical clustering
  
  final_result$model_compare <- model_compare
  final_result$model_labels <- results_list
  final_result$VIres <- vicatmix_labels
  final_result$BayesBinMix <- bayesbin_labels
  final_result
}

save(sim1novarsel, file="sim1novarsel.RData")