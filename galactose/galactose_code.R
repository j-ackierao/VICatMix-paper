#GBM code

galactose <- read.csv("GalactoseData.csv", header=TRUE, row.names=1)

library(tidyverse)
library(mcclust.ext)
source("Variational Mixture Model.R")
set.seed(234)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(25)
galactosedataVI <- foreach(i = 1:25) %dorng% {
  start.time <- Sys.time()
  bigclust <- mixturemodel(galactose, 10, 0.01, 2000, 0.000005)
  end.time <- Sys.time()
  bigclust.timetaken <- end.time - start.time
  start.time2 <- Sys.time()
  littleclust <- mixturemodel(yeastdata, 4, 0.01, 2000, 0.000005)
  end.time2 <- Sys.time()
  littleclust.timetaken <- end.time2 - start.time2
  finalresult <- list()
  #want to record cluster labels
  finalresult$Labels10Clust <- bigclust$model$labels
  finalresult$Labels4Clust <- littleclust$model$labels
  finalresult$Clusters10Clust <- bigclust$Cl[length(bigclust$Cl)]
  finalresult$Clusters4Clust <- littleclust$Cl[length(littleclust$Cl)]
  finalresult$ELBO10Clust <- bigclust$ELBO[length(bigclust$ELBO)]
  finalresult$ELBO4Clust <- littleclust$ELBO[length(littleclust$ELBO)]
  finalresult$Time10Clust <- bigclust.timetaken
  finalresult$Time4Clust <- littleclust.timetaken
  finalresult
}

resultforpsm <- lapply(galactosedataVI, "[[", 1)
p1 <- t(matrix(unlist(resultforpsm), 205, 25))
psm <- comp.psm(p1)
VIcomp <- minVI(psm, method = 'comp',max.k = 10)$cl

resultforpsm <- lapply(galactosedataVI, "[[", 2)
p1 <- t(matrix(unlist(resultforpsm), 205, 25))
psm <- comp.psm(p1)
VIcomp2 <- minVI(psm, method = 'comp',max.k = 4)$cl

galactosedataVI$Avg10Clust <- VIcomp
galactosedataVI$Avg4Clust <- VIcomp2

#Plots

realgalactoselabels <- read.csv("geneNamesClass.csv")
galactose <- read.csv("GalactoseData.csv", header=TRUE, row.names=1)

currentAnnotationRow <- data.frame(
  VICatMixAvg = factor(galactosedataVI$Avg4Clust),
  GOclass = factor(realgalactoselabels$class)
)
rownames(currentAnnotationRow) <- rownames(galactose) 
pheatmap(galactose[sort(galactosedataVI$Avg4Clust, index.return = T)$ix,], 
         cluster_rows = F, show_rownames = F, show_colnames = F, 
         color = colorRampPalette(colors = c("steelblue", "white", "indianred"))(3),
         annotation_row = currentAnnotationRow, cluster_cols = F
)
