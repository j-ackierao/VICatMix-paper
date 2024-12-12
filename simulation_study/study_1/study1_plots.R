#Study 1 plots (all in supplement)
library(tidyverse)

#ELBOARIcorr plot (Section S5.1.1)

ELBOARIcorr <- bind_rows(ELBOARIcorr, .id = "data_number")
colnames(ELBOARIcorr) <- c("data_number", "run", "ELBO", "ARI")

ggplot(ELBOARIcorr, aes(x=ELBO, y=ARI, colour=data_number)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE) + theme_bw() + xlab("log-ELBO")

cor.test(ELBOARIcorr[ELBOARIcorr$data_number == 5,]$ELBO, ELBOARIcorr[ELBOARIcorr$data_number == 5,]$ARI)

#alpha experiments (Section S5.1.2) (change alphaexp2 to alphaexp1 for 10 clusters)

dffinalELBOall <- bind_rows(alphaexp2, .id = "run_number")

#Plots ARI vs alpha
avARI <- dffinalELBOall %>%
  group_by(alpha) %>%
  dplyr::summarize(ARImean = mean(ARI, na.rm=TRUE))

ggplot(data = dffinalELBOall, aes(x=alpha, y=ARI))+
  stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(inherit.aes = TRUE) +
  #geom_line(data = avARI, aes(x = alpha, y = ARImean, group = 1), size = 2, linetype = 'dashed') +
  theme_bw()

#Plots no. of clusters vs alpha
avClust <- dffinalELBOall %>%
  group_by(alpha) %>%
  dplyr::summarize(Clustmean = median(Clusters, na.rm=TRUE))

dffinalELBOall$run_number <- factor(dffinalELBOall$run_number, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

ggplot()+
  geom_jitter(data = dffinalELBOall, aes(x=alpha, y=Clusters, colour=run_number), width = 0.1, height = 0.1) +
  geom_hline(yintercept=10, linetype="dashed", color = "red") +
  #geom_line(data = avClust, aes(x = alpha, y = Clustmean, group = 1), linetype = 'dashed') +
  theme_bw() 

#Changing initial number of clusters (Section S5.1.3)

for (i in 1:10){
  Kexp1[[i]] <- bind_rows(Kexp1[[i]], .id = "run_number")
}
Kexp1 <- bind_rows(Kexp1, .id = "data_number")
avclust <- Kexp1 %>%
  group_by(K) %>%
  dplyr::summarize(clustmean = mean(Clusters, na.rm=TRUE))
avARI <- Kexp1 %>%
  group_by(K) %>%
  dplyr::summarize(ARImean = mean(ARI, na.rm=TRUE))


ggplot() +
  geom_count(data = Kexp1, aes(x=K, y=Clusters), width = 0.05, height = 0.05, shape = 19, alpha = 0.5, color = "darkgreen") +
  theme_bw() + xlab("K_init") + geom_hline(yintercept=4, linetype="dashed", color = "gray") +
  geom_line(data = avclust, aes(x = K, y = clustmean, group = 1)) + geom_point(data = avclust, aes(x = K, y = clustmean, group = 1))

ggplot() +
  geom_jitter(data = Kexp1, aes(x=K, y=ARI, colour=data_number), width = 0, height = 0, shape = 19, alpha = 0.5, color = "darkgreen") +
  theme_bw() + xlab("K_init") + geom_vline(xintercept=4, linetype="dashed", color = "gray") +
  geom_line(data = avARI, aes(x = K, y = ARImean, group = 1)) + geom_point(data = avARI, aes(x = K, y = ARImean, group = 1))


for (i in 1:10){
  Kexp2[[i]] <- bind_rows(Kexp2[[i]], .id = "run_number")
}
Kexp2 <- bind_rows(Kexp2, .id = "data_number")
avclust2 <- Kexp2 %>%
  group_by(K) %>%
  dplyr::summarize(clustmean = mean(Clusters, na.rm=TRUE))
avARI2 <- Kexp2 %>%
  group_by(K) %>%
  dplyr::summarize(ARImean = mean(ARI, na.rm=TRUE))


ggplot() +
  geom_count(data = Kexp2, aes(x=K, y=Clusters), width = 0.05, height = 0.05, shape = 19, alpha = 0.5, color = "darkgreen") +
  theme_bw() + xlab("K_init") + geom_hline(yintercept=10, linetype="dashed", color = "gray") +
  geom_line(data = avclust2, aes(x = K, y = clustmean, group = 1)) + geom_point(data = avclust2, aes(x = K, y = clustmean, group = 1)) + ylim(5, 30)

ggplot() +
  geom_jitter(data = Kexp2, aes(x=K, y=ARI), width = 0, height = 0, shape = 19, alpha = 0.5, color = "darkgreen") +
  theme_bw() + xlab("K_init") + geom_vline(xintercept=10, linetype="dashed", color = "gray") +
  geom_line(data = avARI2, aes(x = K, y = ARImean, group = 1)) + geom_point(data = avARI2, aes(x = K, y = ARImean, group = 1)) + ylim(0.2, 1)



#changing hyperparameter a (Section S5.1.4)

for (i in 1:10){
  cvals_sim[[i]] <- bind_rows(cvals_sim[[i]], .id = "run_number")
}
cvals_sim <- bind_rows(cvals_sim, .id = "data_number")
cvals_sim$c <- as.factor(cvals_sim$c)
colnames(cvals_sim)[3] <- "a"

ggplot(data = cvals_sim, aes(x=a, y=ARI))+
  stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(inherit.aes = TRUE) +
  #geom_line(data = avARI, aes(x = alpha, y = ARImean, group = 1), size = 2, linetype = 'dashed') +
  theme_bw()

compare_means(ARI ~ c, data = cvals_sim[,c(1, 3, 4, 5, 6)], paired = TRUE)

