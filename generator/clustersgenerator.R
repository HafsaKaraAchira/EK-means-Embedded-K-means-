if(!require("clusterGeneration")) { # nolint
    install.packages("clusterGeneration") # nolint
    library("clusterGeneration")
}

library("readr")

K <- 10 # clusters
Dim <- 10 # dimension
# c <- 134218   # cluster points number
# sepvalue <- 0.4   # sep value

for( c in c(1342180) ) {    # 134218,1677725,335545,671090,1342180,2013270,2684360
    for( sepvalue in -9:-9 ) {
        print(c)
        folder <- sprintf("generator/%sD/%sN/SEP%.1f/",Dim,c*10,sepvalue/10.0)       # SEP%.1f
        print(folder)

        minmax_scaler <- function(x) {(x - min(x)) / (max(x) - min(x))}

    #####################################################################################################
    

        tmp1 <- genRandomClust(
                       numClust = K,
                       sepVal = sepvalue/10.0 ,
                       numNonNoisy = Dim,
                       clustszind = 1,
                       clustSizeEq = c, 
                       numReplicate = 1,
                       fileName = paste(folder, "datapoints", sep = "")
                    )

        print("data generated done")

    #####################################################################################################

        data <- read_delim(paste(folder, "datapoints_1.dat", sep = "")," ",col_names= TRUE)
        print("data read done")

        data_norm <- minmax_scaler(data)
        print("data norm done")

        data_norm$cluster <- read_delim(paste(folder, "datapoints_1.mem", sep = "")," ",col_names= FALSE) - 1
        print("data clusters read done") 

        head(data_norm$cluster)

        data_shuffle <- data_norm[sample(1:nrow(data_norm)), ]
        print("data shuffle done")

        write_delim(data.frame(data_shuffle[,c(Dim+1)]),paste(folder, "real_classes.csv", sep = ""),quote = "none",col_names=FALSE,delim = "\t")
        print("data clusters write done")
     
        write_delim(data.frame(data_shuffle[,c(1:Dim)]),paste(folder, "points.csv", sep = ""),quote = "none",col_names=FALSE,delim = "\t")
        print("data norm write done")
    }
}

#################################
#################################

# data_min <- min(data)
# data_max <- max(data)

# word <- readline(prompt="Enter a key (extract the mean vectors before): ")
# print(word)

# original_centres <- read.table(paste(folder,"real_centres.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)
# original_centres_norm <- (original_centres - data_min) / (data_max - data_min)
# write.table(original_centres_norm,paste(folder,"real_centres_norm.csv",sep=""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")

#################################
#################################

# plot(x = data_norm$x1, y = data_norm$x2, xlab = "x1", ylab = "x2", col = rainbow(10)[clusters$V1], pch = 4, main = "clusters")

# clusters <- read.table(paste(folder,"datapoints_1.mem",sep=""), header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# head(clusters)
# points <- read.table("./results/points.csv", header = FALSE,  sep = "\t",  stringsAsFactors = FALSE)
# write.table(data_norm,paste(folder, "points.csv", sep = ""),quote = FALSE,col.names=FALSE,row.names=FALSE, sep = "\t")