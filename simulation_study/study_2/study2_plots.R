#Study 2 plots and further statistical analysis
#Load all relevant RData files from simulations

library(ggpubr)
library(tidyverse)

######################
#SIMULATION 1: n = 1000, p = 100, 10 clusters, no var sel

for (i in 1:10){
  psmtests1000[[i]] <- bind_rows(psmtests1000[[i]], .id = "run_number")
}
psmtests1000 <- bind_rows(psmtests1000, .id = "data_number")
psmtests1000a <- psmtests1000[psmtests1000$run_number != 7,]
psmtests1000a$run_num <- as.factor(psmtests1000a$run_num)

#Adding means into the mix and making a grouped box plot
psmmean1 <- bind_rows(lapply(psmtests1000, "[[", 7), .id = "data_number")
run_nums <- c(5, 10, 15, 20, 25, 30)

psmmean1comp <- list()
for (i in 1:6){
  psmmeantemp <- psmmean1 %>% filter(run_num <= run_nums[i]) %>%
    group_by(data_number) %>%
    dplyr::summarize(Mean_ARI = mean(ARI, na.rm=TRUE), Mean_Clust = mean(Clusters, na.rm = TRUE))
  psmmeantemp <- cbind(psmmeantemp, run_num = run_nums[i])
  psmmean1comp[[i]] <- psmmeantemp
}

psmmean1comp <- bind_rows(psmmean1comp)
psmmean1comp <- cbind(psmmean1comp, model = "Grand Mean", run_number = NA, Time = NA)
psmmean1comp$run_num <- as.factor(psmmean1comp$run_num)
colnames(psmmean1comp)[c(2, 3)] <- c("ARI", "Clusters")
psmmean1comp <- bind_rows(psmtests1000a, psmmean1comp)
psmmean1comp$model[psmmean1comp$model == 'VIavg'] <- 'VoIavg'
psmmean1comp$model[psmmean1comp$model == 'VIcomp'] <- 'VoIcomp'

ggplot(psmmean1comp, aes(x = run_num, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()
ggplot(psmmean1comp, aes(x = run_num, y = Clusters, fill = model)) + geom_hline(yintercept = 10, linetype = 'dashed') + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

wilcox.test(ARI ~ run_num, data = psmmean1comp[psmmean1comp$model == "VoIavg",], paired = TRUE)

pairwise.wilcox.test(x = psmmean1comp[psmmean1comp$model == "VoIcomp",]$ARI, g = psmmean1comp[psmmean1comp$model == "VoIcomp",]$run_num, p.adjust.method = "bonf", paired = TRUE)
pairwise.wilcox.test(x = psmmean1comp[psmmean1comp$run_num == "20",]$ARI, g = psmmean1comp[psmmean1comp$run_num == "20",]$model, p.adjust.method = "bonf", paired = TRUE)

######################
#SIMULATION 2: n = 1000, p = 100, 10 clusters, no var sel, uneven
for (i in 1:10){
  psmtestsuneven[[i]] <- bind_rows(psmtestsuneven[[i]], .id = "run_number")
}
psmtestsunevena <- bind_rows(psmtestsuneven, .id = "data_number")

psmtestsunevena <- psmtestsunevena[psmtestsunevena$run_number != 7,]
psmtestsunevena$run_num <- as.factor(psmtestsunevena$run_num)

psmmean3 <- bind_rows(lapply(psmtestsuneven, "[[", 7), .id = "data_number")
run_nums <- c(5, 10, 15, 20, 25, 30)

psmmean3comp <- list()
for (i in 1:6){
  psmmeantemp <- psmmean3 %>% filter(run_num <= run_nums[i]) %>%
    group_by(data_number) %>%
    dplyr::summarize(Mean_ARI = mean(ARI, na.rm=TRUE), Mean_Clust = mean(Clusters, na.rm = TRUE))
  psmmeantemp <- cbind(psmmeantemp, run_num = run_nums[i])
  psmmean3comp[[i]] <- psmmeantemp
}

psmmean3comp <- bind_rows(psmmean3comp)
psmmean3comp <- cbind(psmmean3comp, model = "Grand Mean", run_number = NA, Time = NA)
psmmean3comp$run_num <- as.factor(psmmean3comp$run_num)
colnames(psmmean3comp)[c(2, 3)] <- c("ARI", "Clusters")
psmmean3comp <- bind_rows(psmtestsunevena, psmmean3comp)
psmmean3comp$model[psmmean3comp$model == 'VIavg'] <- 'VoIavg'
psmmean3comp$model[psmmean3comp$model == 'VIcomp'] <- 'VoIcomp'

ggplot(psmmean3comp, aes(x = run_num, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()
ggplot(psmmean3comp, aes(x = run_num, y = Clusters, fill = model)) + geom_hline(yintercept = 10, linetype = 'dashed') + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()


######################
#SIMULATION 3: n = 2000, p = 100, 20 clusters, no var sel
for (i in 1:10){
  psmtests2000[[i]] <- bind_rows(psmtests2000[[i]], .id = "run_number")
}
psmtests2000a <- bind_rows(psmtests2000, .id = "data_number")
psmtests2000a <- psmtests2000a[psmtests2000a$run_number != 7,]
psmtests2000a$run_num <- as.factor(psmtests2000a$run_num)

psmmean2 <- bind_rows(lapply(psmtests2000, "[[", 7), .id = "data_number")
run_nums <- c(5, 10, 15, 20, 25, 30)

psmmean2comp <- list()
for (i in 1:6){
  psmmeantemp <- psmmean2 %>% filter(run_num <= run_nums[i]) %>%
    group_by(data_number) %>%
    dplyr::summarize(Mean_ARI = mean(ARI, na.rm=TRUE), Mean_Clust = mean(Clusters, na.rm = TRUE))
  psmmeantemp <- cbind(psmmeantemp, run_num = run_nums[i])
  psmmean2comp[[i]] <- psmmeantemp
}

psmmean2comp <- bind_rows(psmmean2comp)
psmmean2comp <- cbind(psmmean2comp, model = "Grand Mean", run_number = NA, Time = NA)
psmmean2comp$run_num <- as.factor(psmmean2comp$run_num)
colnames(psmmean2comp)[c(2, 3)] <- c("ARI", "Clusters")
psmmean2comp <- bind_rows(psmtests2000a, psmmean2comp)
psmmean2comp$model[psmmean2comp$model == 'VIavg'] <- 'VoIavg'
psmmean2comp$model[psmmean2comp$model == 'VIcomp'] <- 'VoIcomp'

ggplot(psmmean2comp, aes(x = run_num, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()
ggplot(psmmean2comp, aes(x = run_num, y = Clusters, fill = model)) + geom_hline(yintercept = 20, linetype = 'dashed') + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

##############################
#Variable selection Simulation 1: 
psmtestsvarsel1a <- lapply(psmtestsvarsel1, head, -2)
for (i in 1:10){
  psmtestsvarsel1a[[i]] <- bind_rows(psmtestsvarsel1a[[i]], .id = "run_number")
}
psmtestsvarsel1a <- bind_rows(psmtestsvarsel1a, .id = "data_number")
psmtestsvarsel1a$run_num <- as.factor(psmtestsvarsel1a$run_num)
psmvarselmean <- bind_rows(lapply(psmtestsvarsel1, "[[", 7), .id = "data_number")

psmvarselmean <- psmvarselmean %>%
  group_by(data_number) %>%
  dplyr::summarize(Mean = mean(ARI, na.rm=TRUE))

psmtestsvarsel1a <- merge(x=psmtestsvarsel1a,y=psmvarselmean, 
                          by="data_number", all.x=TRUE)

psmmeanvar1 <- bind_rows(lapply(psmtestsvarsel1, "[[", 7), .id = "data_number")
run_nums <- c(5, 10, 15, 20, 25, 30)

psmmeanvar1comp <- list()
for (i in 1:6){
  psmmeanvartemp <- psmmeanvar1 %>% filter(run_num <= run_nums[i]) %>%
    group_by(data_number) %>%
    dplyr::summarize(Mean_ARI = mean(ARI, na.rm=TRUE), Mean_Clust = mean(Clusters, na.rm = TRUE))
  psmmeanvartemp <- cbind(psmmeanvartemp, run_num = run_nums[i])
  psmmeanvar1comp[[i]] <- psmmeanvartemp
}

psmmeanvar1comp <- bind_rows(psmmeanvar1comp)
psmmeanvar1comp <- cbind(psmmeanvar1comp, model = "Grand Mean", run_number = NA, Time = NA)
psmmeanvar1comp$run_num <- as.factor(psmmeanvar1comp$run_num)
colnames(psmmeanvar1comp)[c(2, 3)] <- c("ARI", "Clusters")
psmmeanvar1comp <- bind_rows(psmtestsvarsel1a, psmmeanvar1comp)
psmmeanvar1comp$model[psmmeanvar1comp$model == 'VIavg'] <- 'VoIavg'
psmmeanvar1comp$model[psmmeanvar1comp$model == 'VIcomp'] <- 'VoIcomp'

ggplot(psmmeanvar1comp, aes(x = run_num, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()
ggplot(psmmeanvar1comp, aes(x = run_num, y = Clusters, fill = model)) + geom_hline(yintercept = 10, linetype = 'dashed') + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

#############
#Variable selection simulation 2

psmtestsvarsel2a <- lapply(psmtestsvarsel2, head, -2)
for (i in 1:10){
  psmtestsvarsel2a[[i]] <- bind_rows(psmtestsvarsel2a[[i]], .id = "run_number")
}
psmtestsvarsel2a <- bind_rows(psmtestsvarsel2a, .id = "data_number")

psmtestsvarsel2a$run_num <- as.factor(psmtestsvarsel2a$run_num)
psmmeanvar2 <- bind_rows(lapply(psmtestsvarsel2, "[[", 7), .id = "data_number")
run_nums <- c(5, 10, 15, 20, 25, 30)

psmmeanvar2comp <- list()
for (i in 1:6){
  psmmeanvartemp <- psmmeanvar2 %>% filter(run_num <= run_nums[i]) %>%
    group_by(data_number) %>%
    dplyr::summarize(Mean_ARI = mean(ARI, na.rm=TRUE), Mean_Clust = mean(Clusters, na.rm = TRUE))
  psmmeanvartemp <- cbind(psmmeanvartemp, run_num = run_nums[i])
  psmmeanvar2comp[[i]] <- psmmeanvartemp
}

psmmeanvar2comp <- bind_rows(psmmeanvar2comp)
psmmeanvar2comp <- cbind(psmmeanvar2comp, model = "Grand Mean", run_number = NA, Time = NA)
psmmeanvar2comp$run_num <- as.factor(psmmeanvar2comp$run_num)
colnames(psmmeanvar2comp)[c(2, 3)] <- c("ARI", "Clusters")
psmmeanvar2comp <- bind_rows(psmtestsvarsel2a, psmmeanvar2comp)
psmmeanvar2comp$model[psmmeanvar2comp$model == 'VIavg'] <- 'VoIavg'
psmmeanvar2comp$model[psmmeanvar2comp$model == 'VIcomp'] <- 'VoIcomp'

ggplot(psmmeanvar2comp, aes(x = run_num, y = ARI, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()
ggplot(psmmeanvar2comp, aes(x = run_num, y = Clusters, fill = model)) + geom_hline(yintercept = 10, linetype = 'dashed') + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw()

############
#Looking at selected variables 

psmvariables <- lapply(psmtestsvarsel1, "[[", 8)
psmvariables <- lapply(psmvariables, function(x) t(matrix(unlist(x), 100, 30)))

psmthresholds <- list()
run_nums <- c(5, 10, 15, 20, 25, 30)
for (k in 1:6){
  
  psmselect <- list()
  for (i in 1:10){ 
    psmvariables2 <- lapply(psmvariables, function(x) x[1:run_nums[k],])
    psmvariables2[[i]][psmvariables2[[i]] > 0.5] <- 1
    psmvariables2[[i]][psmvariables2[[i]] <= 0.5] <- 0
    result <- vector(mode = 'numeric', length = 100)
    for (j in 1:100){
      result[j] <- sum(psmvariables2[[i]][,j]) / run_nums[k]
    }
    psmselect[[i]] <- result
  }
  
  psmselectthreshold <- data.frame(data_num = rep(c(1:10), 2), threshold = rep(c(0.5, 0.95),each=10), Rel = NA, Irrel = NA, run_num = run_nums[k])
  
  for (i in 1:10){
    psmselectthreshold$Rel[i] <- sum(psmselect[[i]][1:75] >= 0.5)
    psmselectthreshold$Irrel[i] <- sum(psmselect[[i]][76:100] < 0.5)
    psmselectthreshold$Rel[i + 10] <- sum(psmselect[[i]][1:75] >= 0.95)
    psmselectthreshold$Irrel[i + 10] <- sum(psmselect[[i]][76:100] < 0.95)
  }
  
  psmthresholds[[k]] <- psmselectthreshold
  
}

psmthresholds <- bind_rows(psmthresholds)
psmthresholds$threshold <- as.factor(psmthresholds$threshold)
colnames(psmthresholds)[2] <- "model"
psmthresholds$run_num <- as.factor(psmthresholds$run_num)

ggplot(psmthresholds, aes(x = run_num, y = Rel, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw() + geom_hline(yintercept = 75, linetype = 'dotted')
ggplot(psmthresholds, aes(x = run_num, y = Irrel, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw() + geom_hline(yintercept = 25, linetype = 'dotted')
###############
#Simulation 2 for variables selected

psmvariables2 <- lapply(psmtestsvarsel2, "[[", 8)
psmvariables2 <- lapply(psmvariables2, function(x) t(matrix(unlist(x), 100, 30)))

psmthresholds2 <- list()
run_nums <- c(5, 10, 15, 20, 25, 30)
for (k in 1:6){
  
  psmselect2 <- list()
  for (i in 1:10){ 
    psmvariables3 <- lapply(psmvariables2, function(x) x[1:run_nums[k],])
    psmvariables3[[i]][psmvariables3[[i]] > 0.5] <- 1
    psmvariables3[[i]][psmvariables3[[i]] <= 0.5] <- 0
    result <- vector(mode = 'numeric', length = 100)
    for (j in 1:100){
      result[j] <- sum(psmvariables3[[i]][,j]) / run_nums[k]
    }
    psmselect2[[i]] <- result
  }
  
  psmselectthreshold <- data.frame(data_num = rep(c(1:10), 2), threshold = rep(c(0.5, 0.95),each=10), Rel = NA, Irrel = NA, run_num = run_nums[k])
  
  for (i in 1:10){
    psmselectthreshold$Rel[i] <- sum(psmselect2[[i]][1:50] >= 0.5)
    psmselectthreshold$Irrel[i] <- sum(psmselect2[[i]][51:100] < 0.5)
    psmselectthreshold$Rel[i + 10] <- sum(psmselect2[[i]][1:50] >= 0.95)
    psmselectthreshold$Irrel[i + 10] <- sum(psmselect2[[i]][51:100] < 0.95)
  }
  
  psmthresholds2[[k]] <- psmselectthreshold
  
}

psmthresholds2 <- bind_rows(psmthresholds2)
psmthresholds2$threshold <- as.factor(psmthresholds2$threshold)
colnames(psmthresholds2)[2] <- "model"
psmthresholds2$run_num <- as.factor(psmthresholds2$run_num)
psmthresholds2 <- bind_rows(psmthresholds2, psmtestsvarsel2a)

ggplot(psmthresholds2, aes(x = run_num, y = Rel, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw() + geom_hline(yintercept = 50, linetype = 'dotted')
ggplot(psmthresholds2, aes(x = run_num, y = Irrel, fill = model)) + stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.76)) + geom_boxplot(inherit.aes = TRUE, fatten = 1) + theme_bw() + geom_hline(yintercept = 50, linetype = 'dotted')

#Robustness to different initialisations (for supplement)
for (i in 1:10){
  whichrunsinpsm[[i]] <- bind_rows(whichrunsinpsm[[i]][1:10], .id = "run_number")
}
whichrunsinpsm <- bind_rows(whichrunsinpsm, .id = "data_number")

whichrunsinpsm %>%
  group_by(data_number) %>%
  dplyr::summarize(Mean_ARI = mean(ogmean, na.rm=TRUE), Var_ARI = var(ARI, na.rm=TRUE))

###
#F1 scores

psmthresholds$F1 <- (2 * psmthresholds$Rel) / ((2 * psmthresholds$Rel) + (25 - psmthresholds$Irrel) + (75 - psmthresholds$Rel))
psmthresholds2$F1 <- (2 * psmthresholds2$Rel) / ((2 * psmthresholds2$Rel) + (25 - psmthresholds2$Irrel) + (75 - psmthresholds2$Rel))

psmthresholds %>% group_by(run_num, model) %>% filter(run_num == 25 | run_num == 30) %>% dplyr::summarize(Mean_F1 = mean(F1, na.rm=TRUE))



