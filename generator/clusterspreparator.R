if(!require("clusterGeneration")) { 
    install.packages("clusterGeneration")
    library("clusterGeneration")
}

library("readr")

# Elevation / quantitative /meters / Elevation in meters
# Aspect / quantitative / azimuth / Aspect in degrees azimuth
# Slope / quantitative / degrees / Slope in degrees
# Horizontal_Distance_To_Hydrology / quantitative / meters / Horz Dist to nearest surface water features
# Vertical_Distance_To_Hydrology / quantitative / meters / Vert Dist to nearest surface water features
# Horizontal_Distance_To_Roadways / quantitative / meters / Horz Dist to nearest roadway
# Hillshade_9am / quantitative / 0 to 255 index / Hillshade index at 9am, summer solstice
# Hillshade_Noon / quantitative / 0 to 255 index / Hillshade index at noon, summer soltice
# Hillshade_3pm / quantitative / 0 to 255 index / Hillshade index at 3pm, summer solstice
# Horizontal_Distance_To_Fire_Points / quantitative / meters / Horz Dist to nearest wildfire ignition points
# Wilderness_Area (4 binary columns) / qualitative / 0 (absence) or 1 (presence) / Wilderness area designation
# Soil_Type (40 binary columns) / qualitative / 0 (absence) or 1 (presence) / Soil Type designation
# Cover_Type (7 types) / integer / 1 to 7 / Forest Cover Type designation

"Elevation","Aspect","Slope","Horizontal_Distance_To_Hydrology","Vertical_Distance_To_Hydrology","Hillshade_9am","Hillshade_Noon","Hillshade_3pm","Horizontal_Distance_To_Fire_Points","Wilderness_Area_0","Wilderness_Area_1","Wilderness_Area_2","Wilderness_Area_3","Soil_Type_0","Soil_Type_1","Soil_Type_2","Soil_Type_3","Soil_Type_4","Soil_Type_5","Soil_Type_6","Soil_Type_7","Soil_Type_8","Soil_Type_9","Soil_Type_10","Soil_Type_11","Soil_Type_12","Soil_Type_13","Soil_Type_14","Soil_Type_15","Soil_Type_16","Soil_Type_17","Soil_Type_18","Soil_Type_19","Soil_Type_20","Soil_Type_21","Soil_Type_22","Soil_Type_23","Soil_Type_24","Soil_Type_25","Soil_Type_26","Soil_Type_27","Soil_Type_28","Soil_Type_29","Soil_Type_30","Soil_Type_31","Soil_Type_32","Soil_Type_33","Soil_Type_34","Soil_Type_35","Soil_Type_36","Soil_Type_37","Soil_Type_38","Soil_Type_39","Cover_Type"


K <- 6 # clusters
# Dim <- 561 # dimension
# c <- 7350 # 134218   # cluster points number

print(c)
folder <- "generator/real_dataset/covertype/"
print(folder)

minmax_scaler <- function(x) {(x - min(x)) / (max(x) - min(x))}

X_file <- "/home/hafsa/Documents/Paper/covtype/covtype.data"

#"/home/hafsa/Documents/Paper/HAR/X_train.txt"
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