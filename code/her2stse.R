library(Seurat)
library(data.table)
library(STutility)
library(zeallot)
library(openxlsx)

setwd('D:/bioinfo/DLSPseq/ref_information/refcode/her2st-master/her2st-master')
meta_data <- read.xlsx("data/clinical_data/10_HER2+_info.xlsx")
rownames(meta_data) <- meta_data$Sample
samples <- list.files(pattern = ".tsv", path = "data/ST-cnts/", full.names = T)
print(as.data.frame(strsplit(samples, split = "/"))[,3])


#names(samples) <- substr(do.call(rbind, strsplit(samples, split = "/"))[, 10], start = 1, stop = 2)
#imgs <- list.files(path = "data/ST-imgs/", recursive = T, full.names = T, pattern = ".jpg")
#names(imgs) <- do.call(rbind, strsplit(imgs, split = "/"))[, 6]
#ids <- names(samples)
#infoTable <- data.frame(samples, imgs = imgs[ids], ids, patient_id = substr(x = ids, start = 1, stop = 1), stringsAsFactors = FALSE)
#infoTable <- cbind(infoTable, meta_data[infoTable$patient_id, ])
#infoTable[, 8:ncol(infoTable)]
#seu.list <- lapply(unique(infoTable$patient_id), function(s) {
#  InputFromTable(infotable = subset(infoTable, patient_id == s), 
#                 min.gene.spots = 20,
#                 min.spot.feature.count = 300,
#                 platform = "1k")
#}) 
