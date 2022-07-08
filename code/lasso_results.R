# import packages ---------------------------------------------------------

library(glmnet)
library(data.table)
library(doMC)
registerDoMC(cores = 12)
library(dplyr)
library(progress)


# import dataset ----------------------------------------------------------

#inter_val_ST: her2st_st_data.rdata 
#inter_val_feature: her2st_features.rdata
#ex_val_ST: four_ST.rdata
#ex_val_feature: fout_feature.rdata

setwd("D:/bioinfo/DLSPseq/data")

# threshold = 0.5?
sample_ST_zero_filter_0.5 = final_norm_exp %>% .[,colMeans(. != 0) >0.5]


# dataset_generate --------------------------------------------------------

internal_train_test_ratio = 0.8
sample_size = nrow(img_features)
train_size = round(sample_size*internal_train_test_ratio)

set.seed(42)
train_id = sample(1:sample_size,train_size,replace = FALSE)
inter_val_id = setdiff(1:sample_size,train_id)

train_x = as.matrix(img_features[train_id,])

internal_val_x = as.matrix(img_features[inter_val_id,])

ex_val_x = as.matrix(four_feature)

val_gene = intersect(colnames(sample_ST_zero_filter_0.5),colnames(four_ST))


# progress bar ------------------------------------------------------------

pb = progress_bar$new(
  format = '[:bar] :current/:total in :elapsedfull eta: :eta',
  total = length(val_gene), clear = FALSE, width = 80
)


# loop model validation ------------------------------------------------

inter_val_cor_list = list()
ex_val_cor_list = list()

for (i in 1:length(val_gene)) {
  pb$tick()
  
  target_gene = val_gene[i]

  train_y = as.matrix(sample_ST_zero_filter_0.5[,target_gene][train_id])
  
  fit = cv.glmnet(train_x,train_y,family = 'gaussian',alpha = 1,nfolds = 5, parallel = TRUE)
  best = fit$lambda.1se
  
  inter_val_predicted = predict(fit,s = best,newx = internal_val_x)
  inter_val_y = as.matrix(sample_ST_zero_filter_0.5[,target_gene][inter_val_id])
  inter_fit = cor.test(inter_val_predicted,inter_val_y)
  inter_val_cor_list[[i]] = data.frame(gene = target_gene,pvalue = inter_fit$p.value,estimate = inter_fit$estimate)
  
  ex_val_predicted = predict(fit,s = best,newx = ex_val_x)
  ex_val_y = as.matrix(four_ST[,target_gene])
  ex_fit = cor.test(ex_val_predicted,ex_val_y)
  ex_val_cor_list[[i]] = data.frame(gene = target_gene,pvalue = ex_fit$p.value,estimate = ex_fit$estimate)
  }


# pred y cor results list --------------------------------------------------------

inter_val_cor_results =  rbindlist(inter_val_cor_list) %>% na.omit(.)
inter_val_cor_results$pvalue_adj = p.adjust(inter_val_cor_results$pvalue,method = 'BH')

#saveRDS(inter_val_cor_results,file = 'her2st_inter_val_cor_results.rds')

ex_val_cor_results =  rbindlist(ex_val_cor_list) %>% na.omit(.)
ex_val_cor_results$pvalue_adj = p.adjust(ex_val_cor_results$pvalue,method = 'BH')

#saveRDS(ex_val_cor_results,file = 'NBE_ex_val_cor_results.rds')


# val fit gene selection --------------------------------------------------

# p_adj < 0.05? estimate > 0.25?
gene_fit = filter(ex_val_cor_results,pvalue_adj<0.05 & estimate>0.25) 
