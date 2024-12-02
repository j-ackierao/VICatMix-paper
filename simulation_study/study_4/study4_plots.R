library(tidyverse)

#Plots for runtimes study. 
#Results for experiments N is output of 'n_runtimes.R', load in RData output 'n_runtimesnew.RData'
#Results for experiments P is output of 'p_runtimes.R', load in RData output 'p_runtimesnew.RData'

###

#N

for (i in 1:20){
  n_runtimesnew[[i]] <- bind_rows(n_runtimesnew[[i]], .id = "data_number")
}
n_runtimesnew <- bind_rows(n_runtimesnew, .id = "run_number")
n_runtimesnew$n_val <- as.factor(n_runtimesnew$n_val)

ggplot()+ 
  geom_point(data = n_runtimesnew, aes(x=n_val, y=log(Time), color = Model), size = 1.5, shape = 4) +
  theme_bw() + ylab("log(Time (mins))") + xlab("N (number of observations)") #Figure S20 a)

ggplot()+ 
  geom_point(data = n_runtimesnew, aes(x=n_val, y=ARI, color = Model), size = 1.5, shape = 4) +
  theme_bw() + ylab("ARI") + xlab("N (number of observations)") + ylim(0, 1) #Figure S21 a)

###

#P example

for (i in 1:20){
  p_runtimesnew[[i]] <- bind_rows(p_runtimesnew[[i]], .id = "data_number")
}
p_runtimesnew <- bind_rows(p_runtimesnew, .id = "run_number")
p_runtimesnew$p_val <- as.factor(p_runtimesnew$p_val)

ggplot()+ 
  geom_point(data = p_runtimesnew, aes(x=p_val, y=log(Time), color = Model), size = 1.5, shape = 4) +
  theme_bw() + ylab("log(Time (mins))") + xlab("P (number of variables)") #Figure S20 b)

ggplot()+ 
  geom_point(data = p_runtimesnew, aes(x=p_val, y=ARI, color = Model), size = 1.5, shape = 4) +
  theme_bw() + ylab("ARI") + xlab("P (number of variables)") + ylim(-0.02, 1) #Figure S21 b)
