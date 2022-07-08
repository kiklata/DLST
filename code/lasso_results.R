# import packages ---------------------------------------------------------

library(glmnet)
library(data.table)
library(dplyr)
library(progress)


# import dataset ----------------------------------------------------------

#training & inter_val_ST: her2st_st_data.rdata 
#training & inter_val_feature: her2st_features.rdata
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
# internal val gene & external val gene intersect
val_gene = intersect(colnames(sample_ST_zero_filter_0.5),colnames(four_ST))


# progress bar ------------------------------------------------------------

pb = progress_bar$new(
  format = '[:bar] :current/:total in :elapsedfull eta: :eta',
  total = length(val_gene), clear = FALSE, width = 80
)


# loop model validation ------------------------------------------------

inter_val_cor_list_pearson = list()
inter_val_cor_list_spearman = list()
ex_val_cor_list_pearson = list()
ex_val_cor_list_spearman = list()
fit_model_list = list()

for (i in 1:length(val_gene)) {

  pb$tick()

  target_gene = val_gene[i]

  train_y = as.matrix(sample_ST_zero_filter_0.5[,target_gene][train_id])
  # training
  fit = cv.glmnet(train_x,train_y,family = 'gaussian',alpha = 1,nfolds = 5, parallel = F)
  best = fit$lambda.1se
  fit_model_list[[i]] = fit
  fit_model_list[[i]]$gene = target_gene
  
  # internal val
  inter_val_predicted = predict(fit,s = best,newx = internal_val_x)
  inter_val_y = as.matrix(sample_ST_zero_filter_0.5[,target_gene][inter_val_id])
  
  # cor.test(method = pearson? or spearman?)
  inter_fit_pearson = cor.test(inter_val_predicted,inter_val_y,method = c('pearson'))
  inter_val_cor_list_pearson[[i]] = data.frame(gene = target_gene,pvalue = inter_fit_pearson$p.value,estimate = inter_fit_pearson$estimate)
  
  inter_fit_spearman = cor.test(inter_val_predicted,inter_val_y,method = c('spearman'))
  inter_val_cor_list_spearman[[i]] = data.frame(gene = target_gene,pvalue = inter_fit_spearman$p.value,estimate = inter_fit_spearman$estimate)
  
  # external val
  ex_val_predicted = predict(fit,s = best,newx = ex_val_x)
  ex_val_y = as.matrix(four_ST[,target_gene])
  ex_fit_pearson = cor.test(ex_val_predicted,ex_val_y,method = c('pearson'))
  ex_val_cor_list_pearson[[i]] = data.frame(gene = target_gene,pvalue = ex_fit_pearson$p.value,estimate = ex_fit_pearson$estimate)
  
  ex_fit_spearman = cor.test(ex_val_predicted,ex_val_y,method = c('spearman'))
  ex_val_cor_list_spearman[[i]] = data.frame(gene = target_gene,pvalue = ex_fit_spearman$p.value,estimate = ex_fit_spearman$estimate)
  
  }

saveRDS(fit_model_list,file = 'fit_model_list.rds')

# pred y cor results list --------------------------------------------------------
inter_val_cor_results_pearson =  rbindlist(inter_val_cor_list_pearson) %>% na.omit(.)
inter_val_cor_results_pearson$pvalue_adj = p.adjust(inter_val_cor_results_pearson$pvalue,method = 'BH')
saveRDS(inter_val_cor_results_pearson,file = 'her2st_inter_val_cor_results_pearson.rds')

ex_val_cor_results_pearson =  rbindlist(ex_val_cor_list_pearson) %>% na.omit(.)
ex_val_cor_results_pearson$pvalue_adj = p.adjust(ex_val_cor_results_pearson$pvalue,method = 'BH')
saveRDS(ex_val_cor_results_pearson,file = 'NBE_ex_val_cor_results_pearson.rds')

inter_val_cor_results_spearman =  rbindlist(inter_val_cor_list_spearman) %>% na.omit(.)
inter_val_cor_results_spearman$pvalue_adj = p.adjust(inter_val_cor_results_spearman$pvalue,method = 'BH')
saveRDS(inter_val_cor_list_spearman,file = 'her2st_inter_val_cor_results_spearman.rds')

ex_val_cor_results_spearman =  rbindlist(ex_val_cor_list_spearman) %>% na.omit(.)
ex_val_cor_results_spearman$pvalue_adj = p.adjust(ex_val_cor_results_spearman$pvalue,method = 'BH')
saveRDS(ex_val_cor_results_spearman,file = 'NBE_ex_val_cor_results_spearman.rds')

# val fit gene selection --------------------------------------------------

#p_adj < 0.05? estimate > 0.1?
gene_fit_pearson = filter(ex_val_cor_results_pearson,pvalue_adj<0.05 & abs(estimate)>0.1) 

gene_fit_spearman = filter(ex_val_cor_results_spearman,pvalue_adj<0.05 & abs(estimate)>0.1) 


# save selected fit model -------------------------------------------------

selected_gene = gene_fit_pearson$gene

fit_model_list_selected = list()

i = 1  
for (list in fit_model_list) {
  if (list$gene %in% selected_gene) {
    fit_model_list_selected[[i]] = list
    i = i+1
  }}

save(fit_model_list_selected,'fit_model_list_selected.rdata')
