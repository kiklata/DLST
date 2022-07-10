library(magick)
img<- image_read("C:/Users/pz929/Desktop/im/TCGA-A2-A0EY-01Z-00-DX1.png")
img_norm = image_read("C:/Users/pz929/Desktop/im/TCGA-A2-A0EY-01Z-00-DX1.jpg")

library(data.table)

location = svs_wide_d$tiles
pixel = list()
k =1 
for (i in 1:length(location)) {
  pixel[[k]] = data.frame(pixel_x = strsplit(location[i],split = '_')[[1]][1], pixel_y = strsplit(strsplit(location[i],split = '_')[[1]][2],split = '[.]')[[1]][1])
  k=k+1
}

position = rbindlist(pixel)

# scale file
scale_w = 112995/4096
scale_h = 91820/3328  
#TCGA-A2-A0EY-01Z-00-DX1	svs_size	(112995, 91820)	thumb_size	(4096, 3328)

position$scale_x = as.numeric(position$pixel_x)/scale_w
position$scale_y = 3328 - as.numeric(position$pixel_y)/scale_h 

pred_exp = svs_wide_d

position = cbind(position,pred_exp[,2:229])

library(ggplot2)
library(RColorBrewer)

norm = image_ggplot(img_norm)+ labs(title = "Original slide")
  

grey_myplot<- image_ggplot(img)

random_sample_tile = position

mycolors<-brewer.pal(9,"RdBu")
mycolors = rev(mycolors)

library(dplyr)

pp = rbind(random_sample_tile[,1:4],random_sample_tile[,1:4],random_sample_tile[,1:4])

aa = data.frame(exp = random_sample_tile$FASN, gene ='FASN')
bb = data.frame(exp = random_sample_tile$ERBB2, gene ='ERBB2')
cc = data.frame(exp = random_sample_tile$CD68, gene ='CD68')

gg = rbind(aa,bb,cc)

dd = cbind(pp,gg)


p =  grey_myplot +
  geom_point(data = dd,aes(x = scale_x,y = scale_y,color = exp),shape = 16,size =0.91,alpha = 0.7) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20, color = "black"))+
  scale_y_continuous(breaks = NULL)+
  scale_x_continuous(breaks = NULL)+
  scale_color_gradientn(colours = mycolors)+
  theme(legend.title = element_blank())+facet_wrap('gene')+labs(x='',y='')

ggsave('TCGA-A2-A0EY-01Z-00-DX1.jpg',p,width = 32,height = 10,dpi = 300)

ggsave('norm.pdf',norm,width = 10,height = 10,dpi = 72)

p1<- grey_myplot +
  geom_point(data = random_sample_tile,aes(x = scale_x,
                                           y = scale_y,color = FASN),shape = 16,size =2,alpha = 0.6) +
  theme_bw() +
  scale_x_continuous(limits = c(0, as.numeric(image_info(img)[2])), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, as.numeric(image_info(img)[3])), expand = c(0, 0)) +
  scale_color_gradientn(colours = mycolors)+
  theme(legend.title = element_blank())+
  labs(title = "FASN")


p2<- grey_myplot +
  geom_point(data = random_sample_tile,aes(x = scale_x,
                                 y = scale_y,color = CD68),size =2,shape = 16,alpha = 0.6) +
  theme_bw() +
  scale_color_gradientn(colours = mycolors)+
  theme(legend.title = element_blank())+
  labs(title = "CD68")

p3<- grey_myplot +
  geom_point(data = random_sample_tile,aes(x = scale_x,
                                           y = scale_y,color = ERBB2),size =2,shape = 16,alpha = 0.6) +
  theme_bw() +
  scale_x_continuous(limits = c(0, as.numeric(image_info(img)[2])), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, as.numeric(image_info(img)[3])), expand = c(0, 0)) +
  scale_color_gradientn(colours = mycolors)+
  theme(legend.title = element_blank())+
  labs(title = "ERBB2")

library(patchwork)
p4 = p2+p1+p3

ggsave('TCGA-A2-A0EY-01Z-00-DX1.pdf',p4,width = 36,height = 6,dpi = 72)

###
###??Щ?Ҷ??ܺ??? ?????ﶼ?????????˵?
