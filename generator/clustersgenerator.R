if(!require("clusterGeneration")) { # nolint
    install.packages("clusterGeneration") # nolint
    library("clusterGeneration")
}

# library("data.table")
library("readr")

size <- 1000 # dataset size in MB

K <- 10 # clusters

Dim <- 2 # dimension

c <- 20   # cluster points number

sepvalue <- 0.2     # sep value

print(c)

folder <- sprintf("generator/CM13,4M_%sD/test/",Dim,sepvalue)       # SEP%.1f

print(folder)

minmax_scaler <- function(x) {(x - min(x)) / (max(x) - min(x))}

tmp1 <- genRandomClust(
               numClust = K,
               sepVal = sepvalue,
               numNonNoisy = Dim,
               clustszind = 1,
               clustSizeEq = c, 
               numReplicate = 1,
               fileName = paste(folder, "datapoints", sep = "")
            )

print("data generated done")

data <- read_delim(paste(folder, "datapoints_1.dat", sep = "")," ",col_names= TRUE)

# data <- read.delim(paste(folder, "datapoints_1.dat", sep = ""), header=TRUE,sep = " ",as.is=TRUE)

print("data read done")

data_norm <- minmax_scaler(data)

print("data norm done")

write.table(data_norm,paste(folder, "points.csv", sep = ""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")

print("data norm write done")

# clusters <- read.table(paste(folder,"datapoints_1.mem",sep=""), header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# head(clusters)
# points <- read.table("./results/points.csv", header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)

#################################
#################################

# data_min <- min(data)
# data_max <- max(data)

# word <- readline(prompt="Enter a key (extract the mean vectors before): ")
# print(word)

# original_centres <- read.table(paste(folder,"original_centres.csv",sep=""),header = FALSE,sep = " ",stringsAsFactors = FALSE)
# original_centres_norm <- (original_centres - data_min) / (data_max - data_min)
# write.table(original_centres_norm,paste(folder,"original_centres_norm.csv",sep=""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")

#################################
#################################

# plot(x = data_norm$x1, y = data_norm$x2, xlab = "x1", ylab = "x2", col = rainbow(10)[clusters$V1], pch = 4, main = "clusters")
