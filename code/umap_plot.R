library(umap)

data <- umap(`TCGA-A2-A0EY-01Z-00-DX1_pred_ST`[,2:229],n_neighbors=7L,min_dist =0.6,
             method = 'umap-learn',
             metric="manhattan")

library(dplyr)

data$layout %>% 
     as.data.frame() %>% 
     setNames(c("umap1", "umap2")) %>% 
     ggplot(aes(umap1, umap2)) +
     geom_point(size = 2) +
     theme_bw()


library(ggplot2)
library(ggsci)

p <- 
  ggplot(data = as.data.frame(data$layout),aes(V1,V2)) + geom_point() +
  ggtitle(paste0('neighbors_',data$config$n_neighbors,
                 '_min_dist_',data$config$min_dist,
                 '_metric_',data$config$metric)) + scale_color_lancet() +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill='transparent', color='black'),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
        plot.title = element_text(hjust = 0.5, vjust =4,size=20,face = "bold"),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black",face="plain",size=15),
        axis.title = element_text(colour = "black",face="plain",size=15),
        axis.line = element_line(colour = "black",size = 0.5),
        axis.ticks = element_line(colour= "black",size=0.5),
        legend.key = element_rect(fill='transparent', color='transparent'),
        legend.key.size = unit(0.2,"inches"),
        legend.text = element_text(colour = "black", face="plain", size=13))

p
