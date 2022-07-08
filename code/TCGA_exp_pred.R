# import package ----------------------------------------------------------

library(glmnet)
library(data.table)

# import data -------------------------------------------------------------

#gene_fit
#fit_model


# sample ------------------------------------------------------------------

id_to_pred = sample(1:nrow(TCGA_feature),500,replace = FALSE)
tile_to_pred = TCGA_feature[id_to_pred,]

# gene predict ------------------------------------------------------------

# FASN

fit = fit_model_list_selected_pearson_p_adj_0.05_esti_0.1[[112]]

best = fit$lambda.1se

gene_exp_pred = predict(fit,s = best,newx = as.matrix(tile_to_pred))

gene_exp_pred_mean = mean(gene_exp_pred)


# loop --------------------------------------------------------------------



# combine gene matrix ------------------------------------------------------

rbindlist()
t()



