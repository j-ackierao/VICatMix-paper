#DNA methylation
methylationclustA <- read.csv('DNA_Methylation_Cluster_130519.csv', row.names = 1)
methylationclust <- methylationclustA[!is.na(methylationclustA$Cluster),]
rownames(methylationclust) <- NULL

methylationclust$TCGA.ID <- substr(methylationclust$TCGA.ID,1,15)
methduplicates <- which(duplicated(methylationclust$TCGA.ID))

for (d in methduplicates){
  otherdups <- which(methylationclust[,1] == methylationclust[d, 1])
  for (i in 1:length(otherdups)) {
    if (methylationclust[d, 2] != methylationclust[otherdups[i], 2]){
      print(d)
      print(otherdups[i])
      print(methylationclust[d, 2])
      print(methylationclust[otherdups[i], 2])
    }
  }
}
#Disrepency - 2061 and 2060 have the same TCGA ID but different clusters when 12


methylationclustA %>%
  filter(str_detect(TCGA.ID, "TCGA-BH-A1FE"))
#TCGA-BH-A1FE-01A-11D-A13K-05 in Cluster 9, TCGA-BH-A1FE-06A-11D-A212-05 in Cluster 2, TCGA-BH-A1FE-11B-14D-A13T-05 is NA
#Unclear what the extra letters refer to

#methduplicates <- append(methduplicates, 2060)
methylationclust <- methylationclust[-methduplicates,]
colnames(methylationclust)[1] <- "ID"

###########################
#miRNA
mirnaclustA <- read.csv('miRNA.k15.txt', sep = "\t")
mirnaclust <- as.data.frame(cbind(Sample = mirnaclustA$Sample, Cluster = mirnaclustA$Cluster, Disease_code = mirnaclustA$Disease_code))
mirnaclust$Sample <- substr(mirnaclust$Sample,1,15)
mirnaduplicates <- which(duplicated(mirnaclust$Sample)) #Duplicates

for (d in mirnaduplicates){
  otherdups <- which(mirnaclust[,1] == mirnaclust[d, 1])
  for (i in 1:length(otherdups)) {
    if (mirnaclust[d, 2] != mirnaclust[otherdups[i], 2]){
      print(d)
      print(otherdups[i])
      print(mirnaclust[d, 2])
      print(mirnaclust[otherdups[i], 2])
    }
  }
}

#List of duplicates
#1726, 1101
#2224, 1340
#3079, 1422
#3251, 1321
#3617, 1842
#3957, 3296
#3958, 2980

mirnaclustA[c(1726, 1101),]
mirnaclustA[c(2224, 1340),]
mirnaclustA[c(3079, 1422),]
mirnaclustA[c(3251, 1321),]
mirnaclustA[c(3617, 1842),]
mirnaclustA[c(3957, 3296),]
mirnaclustA[c(3958, 2980),]

#Remove all of these for avoidance of doubt

mirnaduplicates <- append(mirnaduplicates, c(1101, 1340, 1422, 1321, 1842, 3296, 2980))
mirnaclust <- mirnaclust[-mirnaduplicates,]
colnames(mirnaclust) <- c("ID", "miRNA", "Disease_code")

############################
#RPPA
rppaclust <- read.csv('PanCan11_RBN_SimpleCluster_20130411.csv')
rppaclust$ID <- gsub(".", "-", rppaclust$ID, fixed = TRUE)
which(duplicated(rppaclust$ID)) #No duplicates

#############################
#mRNA
mrnaclustA <- read.csv('PanCan12.3602-corrected-v3.Subtypes.K16.txt', sep = "\t")
mrnaclust <- mrnaclustA[,c('Sample', 'K16')]
mrnaclust$Sample <- substr(mrnaclust$Sample,1,15)
mrnaduplicates <- which(duplicated(mrnaclust$Sample)) #Duplicates

for (d in mrnaduplicates){
  otherdups <- which(mrnaclust[,1] == mrnaclust[d, 1])
  for (i in 1:length(otherdups)) {
    if (mrnaclust[d, 2] != mrnaclust[otherdups[i], 2]){
      print(d)
      print(otherdups[i])
      print(mrnaclust[d, 2])
      print(mrnaclust[otherdups[i], 2])
    }
  }
}
#Disrepancy in 656, 657
mrnaclustA %>%
  filter(str_detect(Sample, "TCGA-BH-A1FE"))
#TCGA-BH-A1FE-01A-11R-A13Q-07, TCGA-BH-A1FE-06A-11R-A213-07 same one as in DNA methylation

#mrnaduplicates <- append(mrnaduplicates, 656)
mrnaclust <- mrnaclust[-mrnaduplicates,]
colnames(mrnaclust)[1] <- "ID"
#############################
#Somatic copy number
scnaclust <- read.csv('SCNA_Cluster_table.txt', sep = "\t")
scnaclust$ID <- substr(scnaclust$sample,1,15)
which(duplicated(scnaclust$ID)) #No duplicates

##############################
#Put all into one dataframe
#Only keep samples in original COCA analysis 
cocaog <- read.csv('CofC.noMut.K13.Hoadley.20130523.txt', sep = '\t')
cocaog$Samples <- substr(cocaog$Samples,1,15)
samples <- as.data.frame(cocaog$Samples)
colnames(samples) <- 'ID'
which(duplicated(samples)) #No duplicates

samples12 <- samples
samples12$ID <- substr(samples$ID,1,12)

clusterdf12 <- left_join(samples12, rppaclust) #Joins the clusters from RPPA to samples via a left_join in tidyverse

clusterdf <- as.data.frame(cbind(samples$ID, clusterdf12[,2])) #Change the names to 15 characters
colnames(clusterdf) <- c("ID", "RPPA")

clusterdf <- clusterdf %>% left_join(scnaclust[,c("ID", "k8")]) %>% 
  left_join(methylationclust) %>% left_join(mrnaclust) %>% left_join(mirnaclust[,c("ID", "miRNA")])
colnames(clusterdf) <- c("ID", "RPPA", "SCNA", "DNA_meth", "mRNA", "miRNA")

###############################
#Create Matrix of Clusters 

N <- dim(clusterdf)[1] #Number of observations
M <- dim(clusterdf)[2]-1 #Number of datasets

# Number of clusters in each dataset
K <- rep(0, M)
for (i in 1:M) {
  K[i] <- length(unique(na.omit(clusterdf[,i + 1])))
}

#Create matrix
moc <- matrix(0, N, sum(K))
rownames(moc) <- clusterdf$ID
colnames(moc) <- rep(0, sum(K))

#moc[n, j] should be 1 if n is in cluster 

count = 0
for (i in 1:length(K)){
  for (ki in 1:K[i]){
    colnames(moc)[count + ki] <- paste0(colnames(clusterdf)[i+1], "_", unique(na.omit(clusterdf[,i+1]))[ki])
  }
  count <- count + K[i]
  print(count)
}

count = 0
for (i in 1:length(K)){
  for (ki in 1:K[i]){
    #Fill in column by column
    moc[,count + ki] <- clusterdf[,i + 1] == unique(na.omit(clusterdf[,i+1]))[ki]
  }
  count <- count + K[i] #cumulative counting
  print(count)
}


save(moc, file = "MatrixOfClusters.RData")

