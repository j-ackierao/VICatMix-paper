#Study 1 plots (all in supplement)

#ELBOARIcorr plot 

ELBOARIcorr <- bind_rows(ELBOARIcorr, .id = "data_number")
colnames(ELBOARIcorr) <- c("data_number", "run", "ELBO", "ARI")

ggplot(ELBOARIcorr, aes(x=ELBO, y=ARI, colour=data_number)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE) + theme_bw() + xlab("log-ELBO")

cor.test(ELBOARIcorr[ELBOARIcorr$data_number == 5,]$ELBO, ELBOARIcorr[ELBOARIcorr$data_number == 5,]$ARI)

#alpha experiments

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

#changing hyperparameter a 

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
