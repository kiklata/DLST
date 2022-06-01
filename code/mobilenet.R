## MobileNetV2 + 2 FC + 1 classifier
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
  layer_dense(units = 1024, activation = "relu") %>%
  layer_dense(units = 1024, activation = "relu") %>%
  layer_dense(units = 2, activation = "softmax")
model %>% compile(
    loss = "categorical_crossentropy",
    optimizer = "adam",
    metrics = "accuracy")

# Load training data ------------------------------------------------------
image_generator <- image_data_generator( rescale=1/255,validation_split = 0.2)
data_root<- "d:/bioinfo/DLSPseq/ref_information/train_Data/train/"
training_data <- flow_images_from_directory(
  directory = data_root, 
  generator = image_generator,
  target_size = c(224, 224), 
  subset = "training")
validation_data <- flow_images_from_directory(
  directory = data_root, 
  generator = image_generator,
  target_size = c(224, 224), 
  subset = "validation")
model %>% fit(
    training_data, 
    steps_per_epoch = training_data$n/training_data$batch_size,
    validation_data = validation_data,epochs = 100)
save_model_tf(object = model, filepath = "V2_2_normal_epoch16")

# Test model --------------------------------------------------------------
model <- load_model_tf(filepath = "V2_2_normal_epoch16")
image_generator1 <- image_data_generator( rescale=1/255)
test_root<-'d:/bioinfo/DLSPseq/ref_information/train_Data/test/'
test_data <- flow_images_from_directory(
  directory = test_root, 
  generator = image_generator1,
  target_size = c(224, 224), 
  shuffle = FALSE)
steps_per_epoch = test_data$n/test_data$batch_size
test_result <- predict_generator(model, test_data,steps = steps_per_epoch)
evaluate(model,test_data)

