library(ggplot2)
library(ggsci)

data = as.data.frame(data)

p <- 
  ggplot(data = data,aes(V1,V2,color = group_info)) + geom_point() +
  ggtitle(paste0('neighbors_',fe_umap$config$n_neighbors,
                 '_min_dist_',fe_umap$config$min_dist,
                 '_metric_',fe_umap$config$metric)) + scale_color_lancet() +
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
