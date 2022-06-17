library(glmnet)
library(data.table)
library(doMC)
registerDoMC(cores = 8)
library(dplyr)
library(progress)

sample = 'BT23567'
sample_feature = filter(four_feature,substring(rownames(four_feature),1,7)==sample)
sample_ST = filter(four_ST,substring(rownames(four_ST),1,7)==sample)

train_test_ratio = 1
sample_size = nrow(sample_feature)
train_size = round(sample_size*train_test_ratio)

train_x = as.matrix(sample_feature[1:train_size,])

cor_list = list()

sample_ST_zero_filter = t(sample_ST)
sample_ST_zero_filter = sample_ST_zero_filter[rowMeans(sample_ST_zero_filter == 0)<0.8,]
sample_ST_zero_filter = t(sample_ST_zero_filter)

pb = progress_bar$new(
  format = '[:bar] :current/:total in :elapsedfull eta: :eta',
  total = ncol(sample_ST_zero_filter), clear = FALSE, width = 80
)

for (i in 1:ncol(sample_ST_zero_filter)) {

  pb$tick()
  
  target_gene = colnames(sample_ST_zero_filter)[i]
  train_y = as.matrix(sample_ST_zero_filter[,target_gene][1:train_size])
  
  fit = cv.glmnet(train_x,train_y,family ="gaussian",alpha = 1,nfolds = 10, parallel = TRUE)
  best = fit$lambda.1se
  train_predicted = predict(fit,s=best,newx=train_x)
  pred = cor.test(train_predicted,train_y)
  
  cor_list[[i]] = data.frame(gene = target_gene,pvalue = pred$p.value,estimate = pred$estimate)
  
}

#test_x = as.matrix(sample_feature[train_size:sample_size,])
#test_y = as.matrix(sample_ST[,target_gene][train_size:sample_size])

cor_results =  rbindlist(cor_list)
cor_results = na.omit(cor_results)

saveRDS(cor_results,file = paste0(sample,'lasso_cor_results'))