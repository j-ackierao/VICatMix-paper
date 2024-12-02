#Comparing run times with increasing P

library(VICatMix)
library(mclust)

set.seed(2131035)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(20)

p_list <- c(10, 20, 50, 100, 200, 500, 1000)
p_runtimesnew <- foreach(i = 1:20) %dorng% {
  runtimes_p <- list()
  for (j in 1:7){
    generation <- generateSampleDataBin(1000, 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), p_list[j], 0)
    generation2 <- generateSampleDataBin(1000, 10, c(0.1, 0.2, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1), p_list[j] - (p_list[j]/5), p_list[j]/5)
    data <- generation[[1]]
    truelabels <- generation[[2]]
    data2 <- generation2[[1]]
    truelabels2 <- generation2[[2]]
    
    start.time <- Sys.time()
    clustering <- runVICatMix(data, 20, 0.05, tol = 0.00005)
    end.time <- Sys.time()
    timetaken <- difftime(end.time, start.time, unit = "mins")
    
    start.time <- Sys.time()
    clustering2 <- runVICatMixVarSel(data2, 20, 0.05, tol = 0.00005)
    end.time <- Sys.time()
    timetaken2 <- difftime(end.time, start.time, unit = "mins")
    
    modelresult <- data.frame(Model = c("NVS", "VS"), p_val = p_list[j], ELBO = NA, Clusters = NA, ARI = NA, Time = NA)
    modelresult$ELBO[1] <- clustering$ELBO[length(clustering$ELBO)]
    modelresult$Clusters[1] <- clustering$Cl[length(clustering$Cl)]
    modelresult$ARI[1] <- adjustedRandIndex(truelabels$Cluster, clustering$model$labels)
    modelresult$Time[1] <- timetaken
    modelresult$ELBO[2] <- clustering2$ELBO[length(clustering2$ELBO)]
    modelresult$Clusters[2] <- clustering2$Cl[length(clustering2$Cl)]
    modelresult$ARI[2] <- adjustedRandIndex(truelabels2$Cluster, clustering2$model$labels)
    modelresult$Time[2] <- timetaken2
    runtimes_p[[j]] <- modelresult
  }
  runtimes_p
}


save(p_runtimesnew, file="/home/jr951/VariationalMixtures/p_runtimesnew.RData")