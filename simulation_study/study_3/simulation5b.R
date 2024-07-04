#Simulation 2 (variable selection - 50/50, n = 2000 for fun)

library(tidyverse)
library(mclust)
library(BHC)
library(factoextra)
library(stats)
library(mcclust)
library(mcclust.ext)
source("Variational Mixture Model.R")
source("VariationalMixtureModelVarSel.R")
load("~/VariationalMixtures/sim2varsel.RData")
set.seed(1931253)

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
sim2varselBHC <- foreach(i = 1:10) %dorng% {
  generation <- sim2varsel[[i]]$generation
  data <- generation[[1]]
  truelabels <- generation[[2]]
  itemLabels <- rownames(data)
  colnames(data) <- paste0("Covariate",1:100)
  results_list <- list()
  model_compare <- data.frame(model = "BHC", 
                              ARI = NA, Clusters = NA, Time = NA, Rel = NA, Irrel = NA, Rel5 = NA, Irrel5 = NA, RelMean = NA, IrrelMean = NA, RelMean5 = NA, IrrelMean5 = NA)
  #Try every model
  
  #BHC
  start.time <- Sys.time()
  bhcmodel <- bhc(data, itemLabels, verbose=TRUE)
  end.time <- Sys.time()
  bhcclust <- getBHCclusters(bhcmodel, itemLabels)$clusters
  results_list$BHCLabels <- bhcclust
  model_compare$ARI[1] <- adjustedRandIndex(bhcclust, truelabels$Cluster)
  model_compare$Clusters[1] <- length(unique(bhcclust))
  model_compare$Time[1] <- difftime(end.time, start.time, unit = "secs")
  
  results_list$labels <- bhcclust
  results_list$model_compare <- model_compare
  results_list
}

save(sim2varselBHC, file="/home/jr951/VariationalMixtures/sim2varselBHC.RData")