#N
#Replace n_runtimes with n_runtimesvarsel for variable selection

for (i in 1:10){
  n_runtimes[[i]] <- bind_rows(n_runtimes[[i]], .id = "data_number")
}
n_runtimes <- bind_rows(n_runtimes, .id = "run_number")
n_runtimes$n_val <- as.factor(n_runtimes$n_val)
n_runtimes <- n_runtimes[n_runtimes$n_val != 50000,]

ggplot()+ 
  geom_point(data = n_runtimes, aes(x=n_val, y=log(Time)), size = 1.5, shape = 4) +
  theme_bw() + ylab("log(Time (mins))") + xlab("N (number of observations)")

ggplot()+ 
  geom_point(data = n_runtimes, aes(x=n_val, y=ARI), size = 1.5, shape = 4) +
  theme_bw() + ylab("ARI") + xlab("N (number of observations)") + ylim(0, 1)

###

#P example

for (i in 1:10){
  p_runtimes[[i]] <- bind_rows(p_runtimes[[i]], .id = "data_number")
}
p_runtimes <- bind_rows(p_runtimes, .id = "run_number")
p_runtimes$p_val <- as.factor(p_runtimes$p_val)

ggplot()+ 
  geom_point(data = p_runtimes, aes(x=p_val, y=log(Time), color = Model), size = 1.5, shape = 4) +
  theme_bw() + ylab("log(Time (mins))") + xlab("P (number of variables)")