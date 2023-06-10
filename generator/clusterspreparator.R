if(!require("clusterGeneration")) { 
    install.packages("clusterGeneration")
    library("clusterGeneration")
}

library("readr")


K <- 6 # clusters
# Dim <- 561 # dimension
# c <- 7350 # 134218   # cluster points number

print(c)
folder <- "generator/real_dataset/covertype/"
print(folder)

minmax_scaler <- function(x) {(x - min(x)) / (max(x) - min(x))}

X_file <- "/home/hafsa/Documents/Paper/covtype/covtype.data" #"/home/hafsa/Documents/Paper/HAR/X_train.txt"
# Y_file <- "/home/hafsa/Documents/Paper/HAR/y_train.txt"

#####################################################################################################

data <- read_delim(X_file,",",col_names= FALSE)
print("data read done")
head(data)

# data_norm <- minmax_scaler(data[,1:ncol(data)-1])
# print("data norm done")

# data_norm$cluster <- data[,ncol(data)] - 1   # read_delim(Y_file ," ", col_names= FALSE) - 1
# print("data clusters read done")

# head(data_norm$cluster)

# data_shuffle <- data_norm[sample(1:nrow(data_norm)), ]
# print("data shuffle done")

# # [,c(ncol(data_norm)+1)]

# write_delim(data.frame(data_shuffle$cluster),paste(folder, "real_classes.csv", sep = ""),quote = "none",col_names=FALSE,delim = "\t")
# print("data clusters write done")

# write_delim(data.frame(data_shuffle[,c(1:ncol(data_norm)-1)]),paste(folder, "points.csv", sep = ""),quote = "none",col_names=FALSE,delim = "\t")
# print("data norm write done")

# #################################
# #################################

# data_min <- min(data)
# print(data_min)
# data_max <- max(data)
# print(data_max)
# write.table(c(data_min,data_max),paste(folder,"data_min_max.csv",sep=""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")


#################################
#################################

        # # word <- readline(prompt="Enter a key (extract the mean vectors before): ")
        # # print(word)

        # min_max <- read.table(paste(folder,"data_min_max.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)

        # data_min <- min_max[1,1]
        # print(data_min)
        # data_max <- min_max[2,1]
        # print(data_max)

        # original_centres <- read.table(paste(folder,"real_centers.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)
        # original_centres_norm <- (original_centres - data_min) / (data_max - data_min)
        # write.table(original_centres_norm,paste(folder,"real_centers_norm.csv",sep=""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")






#################################
#################################

# plot(x = data_norm$x1, y = data_norm$x2, xlab = "x1", ylab = "x2", col = rainbow(10)[clusters$V1], pch = 4, main = "clusters")

# clusters <- read.table(paste(folder,"datapoints_1.mem",sep=""), header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# head(clusters)
# points <- read.table("./results/points.csv", header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# write.table(data_norm,paste(folder, "points.csv", sep = ""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")
# clusters <- read.table(paste(folder,"datapoints_1.mem",sep=""), header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# head(clusters)
# points <- read.table("./results/points.csv", header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# write.table(data_norm,paste(folder, "points.csv", sep = ""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")