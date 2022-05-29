# H&E Staining to Spatial Transcriptomics using Deep Learning

## Using pre-trained MobileNetV2 model architecture

In R environment, we use interface of tensorflow and keras

```R

setwd() # set workdir

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
data_root<- "/data/train/" # path to data
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
    validation_data = validation_data,epochs = 15)
save_model_tf(object = model, filepath = "V2_2_normal_epoch16")

# Test model --------------------------------------------------------------
model <- load_model_tf(filepath = "V2_2_normal_epoch16")
image_generator1 <- image_data_generator( rescale=1/255)
test_root<-'/data/test/' #path to test
test_data <- flow_images_from_directory(
  directory = test_root, 
  generator = image_generator1,
  target_size = c(224, 224), 
  shuffle = FALSE)
steps_per_epoch = test_data$n/test_data$batch_size
test_result <- predict_generator(model, test_data,steps = steps_per_epoch)
evaluate(model,test_data)
```

## Image normalization

Using adopted code from github

### Training and validation image normalization

We use spatial transcriptomics data (HER2 breast cancer) from her2st<sup>1</sup> as training and internal validation dataset, and ST-net<sup>2</sup> as external validation dataset

### TCGA SVS image normalization

Then we applied our model to TCGA breast cancer H&E image

```python
# coding:utf-8

'''import module'''
from __future__ import division
import os
from pickle import TRUE
#os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import sys
import time
import cv2
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
import openslide 
from resizeimage import resizeimage
import string

print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'started')

args_list = sys.argv

'''imgfile path'''
work_to_be_done = str(args_list[1])
svspath = "/data/TCGA_image/"+work_to_be_done+'/'
resizepath = '/data/resize_img/'+work_to_be_done+'/'
normpath = '/data/norm_img/'+work_to_be_done+'/'
thumbnailpath = '/data/thumbnail/'+work_to_be_done+'/'

os.makedirs(resizepath,exist_ok=True)
os.makedirs(normpath,exist_ok=True)
os.makedirs(thumbnailpath,exist_ok=True)

'''Get magnitude of gradient for given image
ref: https://github.com/yufu2015/PC-CHiP/blob/3145a9e08f90dcdd919fa26c84f4773cedb486ae/inception/preprocess/imgconvert.py'''
def getGradientMagnitude(im):
    im=cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    ddepth = cv2.CV_32F
    dx = cv2.Sobel(im, ddepth, 1, 0)
    dy = cv2.Sobel(im, ddepth, 0, 1)
    dxabs = cv2.convertScaleAbs(dx)
    dyabs = cv2.convertScaleAbs(dy)
    mag = cv2.addWeighted(dxabs, 0.5, dyabs, 0.5, 0)
    return mag

''' normalizeStaining
Source code adapted from: https://github.com/schaugf/HEnorm_python
Original implementation: https://github.com/mitkovetta/staining-normalization'''
def normalizeStaining(img, saveFile=None, Io=240, alpha=1, beta=0.15,sample='sample',img_file='img_file'):             
    HERef = np.array([[0.5626, 0.2159],
                      [0.7201, 0.8012],
                      [0.4062, 0.5581]])
    maxCRef = np.array([1.9705, 1.0308])
    
    # define height and width of image
    h, w, c = img.shape
    # reshape image
    img = img.reshape((-1,3))
    # calculate optical density
    OD = -np.log((img.astype(float)+1)/Io)
    # remove transparent pixels
    ODhat = OD[~np.any(OD<beta, axis=1)]
    # compute eigenvectors
    eigvals, eigvecs = np.linalg.eigh(np.cov(ODhat.T))
    #eigvecs *= -1
    #project on the plane spanned by the eigenvectors corresponding to the two 
    # largest eigenvalues    
    That = ODhat.dot(eigvecs[:,1:3])   
    phi = np.arctan2(That[:,1],That[:,0])
    minPhi = np.percentile(phi, alpha)
    maxPhi = np.percentile(phi, 100-alpha)
    vMin = eigvecs[:,1:3].dot(np.array([(np.cos(minPhi), np.sin(minPhi))]).T)
    vMax = eigvecs[:,1:3].dot(np.array([(np.cos(maxPhi), np.sin(maxPhi))]).T)
    # a heuristic to make the vector corresponding to hematoxylin first and the 
    # one corresponding to eosin second
    if vMin[0] > vMax[0]:
        HE = np.array((vMin[:,0], vMax[:,0])).T
    else:
        HE = np.array((vMax[:,0], vMin[:,0])).T
    # rows correspond to channels (RGB), columns to OD values
    Y = np.reshape(OD, (-1, 3)).T
    # determine concentrations of the individual stains
    C = np.linalg.lstsq(HE,Y, rcond=None)[0]
    # normalize stain concentrations
    maxC = np.array([np.percentile(C[0,:], 99), np.percentile(C[1,:],99)])
    tmp = np.divide(maxC,maxCRef)
    C2 = np.divide(C,tmp[:, np.newaxis])
    # recreate the image using reference mixing matrix
    Inorm = np.multiply(Io, np.exp(-HERef.dot(C2)))
    Inorm[Inorm>255] = 254
    Inorm = np.reshape(Inorm.T, (h, w, 3)).astype(np.uint8)  
    # unmix hematoxylin and eosin
    H = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,0], axis=1).dot(np.expand_dims(C2[0,:], axis=0))))
    H[H>255] = 254
    H = np.reshape(H.T, (h, w, 3)).astype(np.uint8)
    E = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,1], axis=1).dot(np.expand_dims(C2[1,:], axis=0))))
    E[E>255] = 254
    E = np.reshape(E.T, (h, w, 3)).astype(np.uint8)
    if saveFile is not None:
        os.makedirs(normpath+sample,exist_ok=True)
        Image.fromarray(Inorm).save(normpath+sample+'/'+img_file,quality = 100)
    return Inorm, H, E

'''resize TCGA svs to 224X224 pixel jpeg in 40X scope'''
files= os.listdir(svspath)
for filename in files:
    print(filename[0:23])
    img = openslide.open_slide(svspath+str(filename))
    svs_size = img.dimensions
    #generate thumbnail
    thumb = img.get_thumbnail((4096,4096))
    thumb_size = thumb.size
    thumb.save(thumbnailpath + str(filename[0:23]) + '.jpg', 'JPEG', optimize=True, quality=100)
    with open(thumbnailpath+'img_scale.txt','a') as log:
        log.write(str(filename[0:23])+'\t'+'svs_size'+'\t'+str(svs_size)+'\t'+'thumb_size'+'\t'+str(thumb_size)+'\n')
        
    if str(img.properties.values.__self__.get('tiff.ImageDescription')).find('AppMag = 40') != -1:
        sz=220
        seq=400
    else:
        sz=110
        seq=200

    [w, h] = img.dimensions
    os.makedirs(resizepath +str(filename[0:23]),exist_ok=True)
    for x in range(1, w, seq):
        for y in range(1, h, seq):
            img1=img.read_region(location=(x,y), level=0, size=(sz,sz))
            img11=img1.convert("RGB")
            img111=img11.resize((224,224),Image.ANTIALIAS)
            grad=getGradientMagnitude(np.array(img111))
            unique, counts = np.unique(grad, return_counts=True)
            if counts[np.argwhere(unique<=20)].sum() < 224*224*0.6:                
                img111.save(resizepath + str(filename[0:23]) + "/" + str(x) + "_" + str(y) + '.jpg', 'JPEG', optimize=True, quality=100)
                
print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'resized')

'''normalized'''
files= os.listdir(resizepath) 
for filenames in files:
    print(filenames)
    for imgs in os.listdir(resizepath+filenames):
        img =np.array(Image.open(resizepath+filenames+"/"+imgs))
        try:
            normalizeStaining(img = img,
                          saveFile = 1,
                          alpha = 1,
                          beta =0.15,sample=filenames,img_file=imgs)
        except:
            with open(normpath+'error_img.txt','a') as log:
                log.write(filenames+'/'+imgs+'\n')
        
print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'normalized')
```

And process 3X10 svs image in parallel

```shell
# !/bin/bash
python ~/DLST/mycode/preprocessing.py $1 &
python ~/DLST/mycode/preprocessing.py $2 &
python ~/DLST/mycode/preprocessing.py $3
```

save as 'run.sh', then run

```shell
run.sh 20220501 20220502 2220503
```

## Calculated tiles features

```R
# Import Packages ---------------------------------------------------------
library(keras)
library(tfhub)
library(data.table)

# Calculate Feature --------------------------------------------------------------
model = load_model_tf(filepath = 'V2_2_normal_epoch16')
#summary(model)
model_feature = keras_model(inputs = model$input,outputs = get_layer(model, 'global_average_pooling2d_1')$output)

path = '/data/norm_img/20220501' # path to parent dir of pics
sample_list = list.dirs(path, recursive = F)
sample_step = nchar(path) + 2
img_step = nchar(sample_list[1]) + 2
for (i in 1:length(sample_list)) {
  img_list = list.files(sample_list[i], full.names = T)
  print(substring(sample_list[i],sample_step))
  pb = utils::txtProgressBar(max = length(img_list),style = 3)
  img_features = list()
  start_time = Sys.time()
  for (k in 1:length(img_list)) {
    img = image_load(img_list[k], target_size = c(224,224))
    x = image_to_array(img)
    features = predict(model_feature,
                       imagenet_preprocess_input(array_reshape(x, c(1, dim(x)))))
    img_features[[k]] = features
    names(img_features)[[k]] = substring(img_list[k],img_step)
    utils::setTxtProgressBar(pb,k)
  }
  close(pb)
  img_features = as.data.frame(t(rbindlist(list(img_features))))
  save(img_features,file = paste0(sample_list[i],'_train_features.Rdata'))
  end_time = Sys.time()
  print(end_time-start_time)
}
```

run multiple scripts in the background (rstudio jobs)

## File fold example
-DLST
--code
--data
---TCGA_image
---resize_img
---norm_img
---thumbnail
---V2_2_normal_epoch16
--result
--ref_information

## Reference

<sup>1</sup> Spatial deconvolution of HER2-positive breast cancer delineates tumor-associated cell type interactions    
<sup>2</sup> Integrating spatial gene expression and breast tumour morphology via deep learning
