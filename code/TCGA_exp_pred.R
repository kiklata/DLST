# import package ----------------------------------------------------------

library(glmnet)
library(data.table)

# import data -------------------------------------------------------------

#gene_fit
#fit_model


# sample ------------------------------------------------------------------

id_to_pred = sample(1:nrow(TCGA_feature),nrow(TCGA_feature)/20,replace = FALSE)
tile_to_pred = TCGA_feature[id_pred,]

# gene predict ------------------------------------------------------------

best = fit$lambda.1se
gene_exp_pred = predict(fit,s = best,newx = tile_to_pred) 

gene_exp_pred = apply(gene_exp_pred,2,mean)

mean_gene_exp_pred = gene_exp_pred[nrow(gene_exp_pred),]


# loop --------------------------------------------------------------------



# combine gene matrix ------------------------------------------------------

rbindlist()
t()



