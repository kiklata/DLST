library(fitdistrplus)
library(dplyr)

i = 12

print(colnames(four_ST)[i])
dist_test = 'pois'
sample = 'BT23567'

sample_ST = filter(four_ST,substring(rownames(four_ST),1,7)==sample)

sample_ST_zero_filter = t(sample_ST)
sample_ST_zero_filter = sample_ST_zero_filter[rowMeans(sample_ST_zero_filter == 0)<0.2,]
sample_ST_zero_filter = t(sample_ST_zero_filter)

fitd = fitdist(as.numeric(sample_ST_zero_filter[,i]),
               distr = 'norm',method = 'mle')
summary(fitd)
plot(fitd)
