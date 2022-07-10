rm(list = ls())

library(glmnet)
library(data.table)
library(reshape2)  

setwd('D:/bioinfo/DLSPseq/data/features')

#slide_list = list.files(getwd())

#slide_exp = list()


# single slide predict

svs = 'TCGA-A2-A0EY-01Z-00-DX1'

tile_to_pred = img_features

fit_length = length(fit_model_list_selected_pearson_p_adj_0.05_esti_0.1)

svs_exp_pred = list()
for (i in 1:fit_length) {

  fit = fit_model_list_selected_pearson_p_adj_0.05_esti_0.1[[i]]
  gene = fit$gene
  best = fit$lambda.1se
  exp_pred = predict(fit,s = best,newx = as.matrix(tile_to_pred))
  
  svs_exp_pred[[i]] = data.frame(genes = gene, exp = exp_pred,tiles = rownames(exp_pred))
  
}

svs_long_m = rbindlist(svs_exp_pred)
svs_wide_d<-dcast(svs_long_m, tiles~genes,value.var = 's1')

saveRDS(svs_wide_d,file = paste0(svs,"_pred_ST.rds"))

# LOOP
for (k in 1:length(slide_list)) {
  
  slide = slide_list[k]
  
  slide_name = substring(slide,1,23)
  
  load(paste0(slide))

  if (nrow(img_features)>500) {
    # if all feature calculated  
    id_to_pred = sample(1:nrow(img_features),500,replace = FALSE)
    tile_to_pred = img_features[id_to_pred,]
  }
  else{
    # if 500 feature calculated
    tile_to_pred = img_features
  }
    

  fit_length = length(fit_model_list_selected_pearson_p_adj_0.05_esti_0.1)
  
  mean_exp_pred = list()
  for (i in 1:fit_length) {
    fit = fit_model_list_selected_pearson_p_adj_0.05_esti_0.1[[i]]
    gene = fit$gene
    best = fit$lambda.1se
    exp_pred = predict(fit,s = best,newx = as.matrix(tile_to_pred))
    mean_exp_pred[[i]] = data.frame(genes = gene, mean_exp = mean(exp_pred),slides = slide_name)
  }
  
  slide_exp[[k]] = rbindlist(mean_exp_pred)
  
  }



# combine gene matrix ------------------------------------------------------

data_long_m = rbindlist(slide_exp)

data_wide_d<-dcast(data_long_m, genes~slides,value.var = 'mean_exp')


# TCGA cor test -----------------------------------------------------------

pearson_cor_res_list = list()
spearman_cor_res_list = list()

for (i in 1:nrow(myy)) {
  
  gene_name = data$genes[i]
  x = as.numeric(data[i,][2:ncol(data)])
  y = as.numeric(myy[i,][2:ncol(myy)])
  
  res_pearson = cor.test(x,y,method = c('pearson'))
  res_spearman = cor.test(x,y,method = c('spearman'))
  
  pearson_cor_res_list[[i]] = data.frame(gene = gene_name,cor = res_pearson$estimate,p = res_pearson$p.value)
  spearman_cor_res_list[[i]] = data.frame(gene = gene_name,cor = res_spearman$estimate,p = res_spearman$p.value)
  
}

