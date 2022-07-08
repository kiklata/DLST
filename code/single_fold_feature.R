# Import Packages ---------------------------------------------------------
library(keras)
library(tfhub)
library(data.table)
library(progress)
tensorflow::set_random_seed(42,disable_gpu = F)

# Calculate Feature --------------------------------------------------------------
model = load_model_tf(filepath = 'D:/bioinfo/DLSPseq/data/V2_2_normal_epoch16')
#summary(model)
model_feature = keras_model(inputs = model$input,outputs = get_layer(model, 'global_average_pooling2d_1')$output)

path = 'D:/bioinfo/DLSPseq/data/he2st_st_data/her2st_normed/tile'
img_list = list.files(path, full.names = T)

pb = progress_bar$new(
    format = '[:bar] :current/:total in :elapsedfull eta: :eta',
    total = length(img_list), clear = FALSE, width = 80)
  
img_features = list()
for (k in 1:length(img_list)) {
    img = image_load(img_list[k], target_size = c(224,224))
    x = image_to_array(img)
    features = predict(model_feature,
                       imagenet_preprocess_input(array_reshape(x, c(1, dim(x)))))
    img_features[[k]] = features
    names(img_features)[[k]] = substring(img_list[k],58)
    
    pb$tick()
    
  }
img_features = as.data.frame(t(rbindlist(list(img_features))))
save(img_features,file = paste0('D:/bioinfo/DLSPseq/her2st_features.Rdata'))


