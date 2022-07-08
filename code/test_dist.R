library(fitdistrplus)
library(dplyr)

i = 12

print(colnames(four_ST)[i])

sample_ST_zero_filter_0.5 = final_norm_exp %>% .[,colMeans(. != 0) >0.5]

dim(sample_ST_zero_filter_0.5)

fitd = fitdist(as.numeric(sample_ST_zero_filter_0.5[,1]),
               distr = 'norm',method = 'mle')
summary(fitd)
plot(fitd)
