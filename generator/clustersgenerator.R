if(!require("clusterGeneration")) { # nolint
    install.packages("clusterGeneration") # nolint
    library("clusterGeneration")
}

size <- 1000 # dataset size in MB

K <- 10 # clusters

D <- 10 # dimension

ligne_avg_bytes <- D * 17 + (D-1)

folder <- "generator/CM5,8M_1000MO_SEP0,3/" 

# "generator/CM4_5,8M_1000MO/"

c <- 580000 #( signif(size / ligne_avg_bytes , digits = 2) * 1000000 ) / K   # 580000

print(c)

minmax_scaler <- function(x) {(x - min(x)) / (max(x) - min(x))}

# tmp1 <- genRandomClust(
#                numClust = K,
#                sepVal = 0.3,
#                numNonNoisy = D,
#                clustszind = 1,
#                clustSizeEq = c, 
#                numReplicate = 1,
#                fileName = paste(folder, "datapoints", sep = "")
#                 )

data <- read.delim(paste(folder,"datapoints_1.dat",sep=""), header=TRUE,sep = " ",as.is=TRUE)

# data_norm <- minmax_scaler(data)

# write.table(data_norm,paste(folder,"points.csv",sep=""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")


# clusters <- read.table(paste(folder,"datapoints_1.mem",sep=""), header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# head(clusters)
# points <- read.table("./results/points.csv", header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)

#################################
#################################

data_min <- min(data)
data_max <- max(data)

word <- readline(prompt="Enter a key (extract the mean vectors before): ")
print(word)

original_centres <- read.table(paste(folder,"original_centres.csv",sep=""),header = FALSE,sep = " ",stringsAsFactors = FALSE)
original_centres_norm <- (original_centres - data_min) / (data_max - data_min)
write.table(original_centres_norm,paste(folder,"original_centres_norm.csv",sep=""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")

#################################
#################################

# plot(x = data_norm$x1, y = data_norm$x2, xlab = "x1", ylab = "x2", col = rainbow(10)[clusters$V1], pch = 4, main = "clusters")
