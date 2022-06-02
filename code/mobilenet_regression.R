## multi-output regression model

# Import packages ---------------------------------------------------------
library(keras)
library(tfhub)

# Initialize model --------------------------------------------------------
image_shape <- c(224L, 224L, 3L)
basemodel<-application_mobilenet_v2(weights = 'imagenet',include_top = F,input_shape = image_shape) 
freeze_weights(basemodel)
model <- keras_model_sequential() %>%
  basemodel %>%
  layer_global_average_pooling_2d()%>%
  layer_dropout(rate = 0.6)%>%
#  layer_dense(units = 1024, activation = "relu") %>%
#  layer_dense(units = 1024, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")
model %>% compile(
  loss = "mse",
  optimizer = "adam",
  metrics = list("mse","mae"))
# Load training data ------------------------------------------------------
image_generator <- image_data_generator(rescale=1/255,validation_split = 0.2)
data_root<- "d:/bioinfo/DLSPseq/data/BC_NG/normed_tiles/"
training_data <- flow_images_from_dataframe(
  ST,
  directory = data_root, 
  x_col = 'ID',
  y_col = 'SDC1',
  generator = image_generator,
  class_mode = 'other',
  target_size = c(224, 224), 
  subset = "training")
validation_data <- flow_images_from_dataframe(
  ST,
  directory = data_root, 
  x_col = 'ID',
  y_col = 'SDC1',
  generator = image_generator,
  class_mode = 'other',
  target_size = c(224, 224), 
  subset = "validation")
model %>% fit(
  training_data, 
  steps_per_epoch = training_data$n/training_data$batch_size,
  validation_data = validation_data,epochs = 20)
save_model_tf(object = model, filepath = "V2_2_normal_epoch16")

# Test model --------------------------------------------------------------
model <- load_model_tf(filepath = "V2_2_normal_epoch16")
image_generator1 <- image_data_generator( rescale=1/255)
test_root<-'d:/bioinfo/DLSPseq/data/BC_NG/normed_tiles/'
test_data <- flow_images_from_dataframe(
  ST,
  directory = data_root, 
  x_col = 'ID',
  y_col = 'SDC1',
  generator = image_generator,
  class_mode = 'other',
  target_size = c(224, 224), 
  subset = "validation")
steps_per_epoch = test_data$n/test_data$batch_size
test_result <- predict(model, test_data,steps = steps_per_epoch)

evaluate(model,test_data)

