#Running VICatMix on pan-cancer data

load("MatrixOfClusters.RData")
source("Variational Mixture Model.R")
library(tidyverse)
library(mclust)
library(mcclust)
library(mcclust.ext)

moc[is.na(moc)] <- 0

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(25)
cocaclustResultAv <- foreach(i = 1:25) %dorng% {
  start.time <- Sys.time()
  cocaclust <- mixturemodel(moc, 15, 0.05, 1500, 0.000005)
  end.time <- Sys.time()
  timetaken <- end.time - start.time
  finalresult <- list()
  #want to record cluster labels
  finalresult$Labels <- cocaclust$model$labels
  finalresult$Clusters <- cocaclust$Cl[length(cocaclust$Cl)]
  finalresult$ELBO <- cocaclust$ELBO[length(cocaclust$ELBO)]
  finalresult$Time <- timetaken
  finalresult
}

labels <- lapply(cocaclustResultAv[1:25], "[[", 1)
p1 <- t(matrix(unlist(labels), 3527, 25))
psm <- comp.psm(p1)
VIavg <- minVI(psm, method = 'comp',max.k = 15)$cl

cocaclustResultAv$AV <- VIavg

###
#Creating plots 

#Heatmap - all samples
cocaav15 <- cocaclustResultAv$AV

#Tissues of origin
cocaog <- read.csv('CofC.noMut.K13.Hoadley.20130523.txt', sep = '\t')

currentAnnotationRow <- data.frame(
  Cluster = factor(cocaav15),
  Tissue = factor(cocaog$Tissue)
)

annotationColor <- list(Tissue  = c("BLCA" = "#f3c300", "BRCA" = "#875692", "COAD" = "#f38400", "GBM" = "#a1caf1", "HNSC" = "#be0032", "KIRC" = "#c2b280",
                                    "LAML" = "#848482", "LUAD" = "#008856", "LUSC" = "#e68fac", "OV" = "#0067a5", "READ" = "#f99379", "UCEC" = "#604e97"),
                        Cluster = c("1" = "#A6CEE3", "2" = "#1F78B4", "3" = "#B2DF8A", "4" = "#33A02C", "5" = "#FB9A99", "6" = "#E31A1C", "7" = "#FDBF6F", 
                                    "8" = "#FF7F00", "9" = "#CAB2D6", "10" = "#6A3D9A", "11" = "#FFFF99", "12" = "#B15928", "13" = "#1ff8ff", "14" = "#1B9E77", "15" = "#D95F02")
)

rownames(currentAnnotationRow) <- rownames(moc)

pheatmap(moc[order(cocaog$Tissue, cocaav15),], 
         show_rownames = F, annotation_row = currentAnnotationRow, 
         annotation_colors= annotationColor,
         color = colorRampPalette(colors = c("white", "black"))(2),
         cluster_rows = F, cluster_cols = F, fontsize = 8)

#Percentage heat map

COCAPercentHeat15 <- matrix(0, nrow = 12, ncol = 15)
tissuesandclusters <- as.data.frame(cbind(cocaog$Tissue, cocaav15))

for (i in 1:12){
  for (j in 1:15){
    COCAPercentHeat15[i, j] <- sum(tissuesandclusters$V1 == tissues[i] & tissuesandclusters$cocaav15 == j) / sum(cocaog$Tissue == tissues[i])
  }
}

rownames(COCAPercentHeat15) <- levels(factor(cocaog$Tissue))
colnames(COCAPercentHeat15) <- paste0(seq(ncol(COCAPercentHeat15)))

pheatmap(COCAPercentHeat15, display_numbers = F, 
         cluster_rows = F, cluster_cols = F, fontsize_number = 10,
         color = colorRampPalette(colors = c("white", "black"))(198),
         show_rownames = T, show_colnames = T, angle_col = 0)


#BRCA heat map

COCAbestlabelsdf15 <- as.data.frame(cbind(cocaog, cocaav15))
PAM50 <- read.csv("TCGAPAM50_2018.csv", header = TRUE)
PAM50 <- PAM50 %>% select(Sample.ID, BRCA_Subtype_PAM50)
BRCAcoca15 <- COCAbestlabelsdf15 %>% filter(Tissue == "BRCA")
BRCAcoca15 <- BRCAcoca15 %>% 
  mutate(Sample.ID = str_extract(Samples, "^.{1,12}"))

BRCAcoca15 <- left_join(BRCAcoca15, PAM50, by = "Sample.ID")

BRCAlist <- which(cocaog$Tissue == "BRCA")

BRCAmoc <- moc[BRCAlist,]

currentAnnotationRow <- data.frame(
  PAM50 = factor(BRCAcoca15$BRCA_Subtype_PAM50),
  Cluster = factor(BRCAcoca15$cocaav15)
)
rownames(currentAnnotationRow) <- rownames(BRCAmoc)

pheatmap(BRCAmoc[order(BRCAcoca15$cocaav15, BRCAcoca15$BRCA_Subtype_PAM50),], fontsize = 7,
         show_rownames = F, annotation_row = currentAnnotationRow, 
         #annotation_colors= annotationColor,
         color = colorRampPalette(colors = c("white", "black"))(2),
         cluster_rows = F, cluster_cols = F, name = "score")

freqtableBRCA <- BRCAcoca15 %>% 
  count(BRCA_Subtype_PAM50, cocaav15) %>% 
  pivot_wider(names_from = cocaav15, values_from = n, values_fill = list(n = 0), names_sort = TRUE)





