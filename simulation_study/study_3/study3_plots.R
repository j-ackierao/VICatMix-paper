#Simulations and comparisons - plots
#Replace 'sim1novarsel_new' with relevant data eg. 'sim2novarsel_new'
library(dplyr)
library(ggplot2)

#No variable selection example
simsplot <- lapply(sim1novarsel_new, "[[", 1) 
simsplot <- bind_rows(simsplot, .id = "data_number")
simsplot <- simsplot[simsplot$model != "FlexMixBIC",] #All the same for BIC and ICL
#Comparison of ARI
ggplot(simsplot, aes(x = model, y = ARI)) + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(inherit.aes = TRUE, fatten = 1, fill = "lavender") + theme_bw() + ylim(0,1)
#Comparison of number of clusters
ggplot(simsplot, aes(x = model, y = Clusters)) +geom_hline(yintercept=10, linetype='dashed', col = 'red') + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(inherit.aes = TRUE, fatten = 1, fill = "lightsalmon") + theme_bw()

#Variable selection example
simsvsplot <- lapply(sim1varsel_new, "[[", 1)
simsvsplot <- bind_rows(simsvsplot, .id = "data_number")
simsvsplot <- simsvsplot[simsvsplot$model != "FlexMixBIC",]
simsvsplot$Rel <- coalesce(simsvsplot$Rel5, simsvsplot$Rel)
simsvsplot$Irrel <- coalesce(simsvsplot$Irrel5, simsvsplot$Irrel)
#Comparison of ARI
ggplot(simsvsplot, aes(x = model, y = ARI)) + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(inherit.aes = TRUE, fatten = 1, fill = "lavender") + theme_bw() + ylim(0,1) + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
#Comparison of number of clusters
ggplot(simsvsplot, aes(x = model, y = Clusters)) +geom_hline(yintercept=10, linetype='dashed', col = 'red') + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(inherit.aes = TRUE, fatten = 1, fill = "lightsalmon") + theme_bw() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#Times
simstimes <- simsplot %>%
  group_by(model) %>%
  dplyr::reframe(x = quantile(Time, c(0.25, 0.5, 0.75), na.rm = TRUE), q = c(0.25, 0.5, 0.75))

simsvstimes <- simsvsplot %>%
  group_by(model) %>%
  dplyr::reframe(x = quantile(Time, c(0.25, 0.5, 0.75), na.rm = TRUE), q = c(0.25, 0.5, 0.75))

#F1 scores
simsvsplot$F1 <- (2 * simsvsplot$Rel) / ((2 * simsvsplot$Rel) + (50 - simsvsplot$Irrel) + (50 - simsvsplot$Rel))
simsvsplot %>% group_by(model) %>% dplyr::summarize(Mean_F1 = mean(F1, na.rm=TRUE))

psmthresholds2$F1 <- (2 * psmthresholds2$Rel) / ((2 * psmthresholds2$Rel) + (25 - psmthresholds2$Irrel) + (75 - psmthresholds2$Rel))


