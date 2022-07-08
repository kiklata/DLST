# coding:utf-8

'''import module'''
from __future__ import division
import os
#os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import sys
import time
import cv2
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from alive_progress import alive_bar


args_list = sys.argv

'''imgfile path'''
work_to_be_done = str(args_list[1])
resizepath = 'd:/bioinfo/DLSPseq/data/resize_img/'+work_to_be_done+'/'
normpath = 'd:/bioinfo/DLSPseq/data/norm_img/'+work_to_be_done+'/'

os.makedirs(normpath,exist_ok=True)

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
    #H = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,0], axis=1).dot(np.expand_dims(C2[0,:], axis=0))))
    #H[H>255] = 254
    #H = np.reshape(H.T, (h, w, 3)).astype(np.uint8)
    #E = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,1], axis=1).dot(np.expand_dims(C2[1,:], axis=0))))
    #E[E>255] = 254
    #E = np.reshape(E.T, (h, w, 3)).astype(np.uint8)
    if saveFile is not None:
        os.makedirs(normpath+sample,exist_ok=True)
        Image.fromarray(Inorm).save(normpath+sample+'/'+img_file,quality = 100)
    return Inorm #, H, E

'''normalized'''
def norm_func():
    with alive_bar(len(os.listdir(resizepath))) as bar:
        for imgs in os.listdir(resizepath):
            bar()
            img =np.array(Image.open(resizepath+"/"+imgs))
            try:
                normalizeStaining(img = img,
                                saveFile = 1,
                                alpha = 1,
                                beta =0.15,sample='',img_file=imgs)
            except:
                with open(normpath+'error_img.txt','a') as log:
                    log.write(imgs+'\n')
            
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),work_to_be_done,' normalized')

norm_func()