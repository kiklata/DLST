library(glmnet)
library(data.table)
library(doMC)
registerDoMC(cores = 8)
library(dplyr)
library(progress)
# library(mpath)
library(ggplot2)

sample = 'BT23567'
sample_feature = filter(four_feature,substring(rownames(four_feature),1,7)==sample)
sample_ST = filter(four_ST,substring(rownames(four_ST),1,7)==sample)

train_test_ratio = 1
sample_size = nrow(sample_feature)
train_size = round(sample_size*train_test_ratio)

train_x = as.matrix(sample_feature[1:train_size,])

sample_ST_zero_filter = t(sample_ST)
sample_ST_zero_filter = sample_ST_zero_filter[rowMeans(sample_ST_zero_filter == 0)<0.9,]
sample_ST_zero_filter = t(sample_ST_zero_filter)


exp = list()
for (i in 1:20) {
  target_gene = colnames(sample_ST_zero_filter)[i]
  train_y = as.matrix(sample_ST_zero_filter[,target_gene][1:train_size])
  exp[[i]] = data.frame(img = rownames(train_x),gene = target_gene,exps = train_y)
  #p = ggplot(as.data.frame(train_y))+geom_histogram(aes(train_y),bins = 10)+labs(title = target_gene)
  
}
exp = rbindlist(exp)
plotlist = list()
genes = unique(exp$gene)
for (g in genes) {
  p = ggplot(exp %>% filter(.,gene == g))+geom_histogram(aes(exps),bins = 10)+labs(title = g)
  plotlist[[g]] = p
}
library(cowplot)
all = plot_grid(plotlist =  plotlist,ncol = 5 ,align = "hv")
hei = length(plotlist)/5*5
ggsave(all,
       filename = 'all.jpg',
       limitsize = FALSE,
       width = 20,
       height = hei)

cor_list = list()
pb = progress_bar$new(
  format = '[:bar] :current/:total in :elapsedfull eta: :eta',
  total = 20, clear = FALSE, width = 80
)
for (i in 1:20) {

  pb$tick()
  target_gene = colnames(sample_ST_zero_filter)[i]
  train_y = as.matrix(sample_ST_zero_filter[,target_gene][1:train_size])

  fit = cv.glmnet(train_x,train_y,family = 'poisson',alpha = 1,nfolds = 10, parallel = TRUE)
#  fit = cv.glmregNB(train_y~train_x)
  
  best = fit$lambda.1se
  train_predicted = predict(fit,s=best,newx=train_x)
  pred = cor.test(train_predicted,train_y)
  
  cor_list[[i]] = data.frame(gene = target_gene,pvalue = pred$p.value,estimate = pred$estimate)

  }

# test_x = as.matrix(sample_feature[train_size:sample_size,])
# test_y = as.matrix(sample_ST[,target_gene][train_size:sample_size])

cor_results =  rbindlist(cor_list)
cor_results = na.omit(cor_results)
cor_results$pvalue_adj = p.adjust(cor_results$pvalue,method = 'BH')

saveRDS(cor_results,file = paste0(sample,'_lasso_cor_results.rds'))
