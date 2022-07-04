# H&E Staining to Spatial Transcriptomics using Deep Learning

## Using pre-trained MobileNetV2 model architecture

In R environment, we use interface of tensorflow and keras

```R

setwd('DLST') # set workdir

source('code/mobilenet_classification')
```

## Image normalization

Using adopted code from github

### Training and validation image normalization

We use spatial transcriptomics data (HER2 breast cancer) from her2st<sup>1</sup> as training and internal validation dataset, and ST-net<sup>2</sup> as external validation dataset

### TCGA SVS image normalization

Then we applied our model to TCGA breast cancer H&E image

```python
code/preprocessing.py 20220501
```

## Calculate tiles features

```R
source('code/feature_calc')
```

## File folder example

```
-DLST 
--code  
---feature_calculate.R  
---mobilenet.R  
---preprocessing.py 
---run.sh 
--data  
---TCGA_image # each fold contains 10 svs image
----20220501 
----20220502
----20220503
---resize_img 
---norm_img 
---thumbnail  
---V2_2_normal_epoch16  
--result  
--ref_information 
--readme.md
```

## normed TCGA image tiles

[Tiles](https://pan.baidu.com/s/1jf_9ckPXsCGwZGybsVfm8A?pwd=zdh5)

## Reference

<sup>1</sup> Spatial deconvolution of HER2-positive breast cancer delineates tumor-associated cell type interactions    
<sup>2</sup> Integrating spatial gene expression and breast tumour morphology via deep learning
