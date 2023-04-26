library(remotes) # nolint

data_norm <- read.table('./results/points_0,56M_100MO.csv', header = FALSE,  sep = '\t',  stringsAsFactors = FALSE) # nolint

knitr::opts_chunk$set(
  dpi = 150,  # nolint
  fig.width = 6, 
  fig.height = 4,
  message = FALSE
  # warning = FALSE
)

data_norm.pca <- prcomp(data_norm[,(1:10)], center = TRUE,scale. = TRUE) # nolint

summary(data_norm.pca)

library(devtools)
library(ggbiplot)

#summary(data_norm.pca)
#str(data_norm.pca)

#biplot(data_norm.pca)
#ggbiplot(data_norm.pca)

data_norm.clusters <- read.table('./generator/datapoints_1.mem', header = FALSE,  sep = '\t',  stringsAsFactors = FALSE)

ggbiplot::ggbiplot(data_norm.pca,obs.scale = 1,var.scale = 1,ellipse = TRUE,circle=TRUE,groups = clusters)
        #var.axes = FALSE,
        #labels = rownames(mtcars),    
#)

#+ scale_colour_manual(name = "Origin",
#                      values = rainbow(10)
#  ) +
#  ggtitle("PCA of clusters dataset") +
#  theme_minimal() +
#  theme(legend.position = "bottom")


