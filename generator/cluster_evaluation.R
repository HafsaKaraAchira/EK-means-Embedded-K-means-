
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) # alternatively, this also loads %>%
library(mclust)
library(clue)

library(pdist)
library("sos")
library("matrixStats")

#############################################
#             clusters centers
#############################################

args <- commandArgs()
print(args)

SEP <- args[6]
N <- args[7]
M <- args[8]
K <- args[9]
DIM <- args[10]
DMAX <- args[11]

########################################
########################################


report <- sprintf("SKIP_CHUNKS/reports/%ldN_%ldM_%ldD_%ldK_%dL/",N,M,DIM,K,DMAX)

folder <- sprintf("generator/%sD/%sN/SEP%s/",DIM,N,SEP) 

original_centres <- read.table(paste(folder,"real_centers_norm.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)
clusters_centres <- read.table(paste(result,"result_centers.csv"), header = FALSE,  sep = "\t", stringsAsFactors = FALSE)

########################################
########################################
########################################
########################################
########################################

#get the real and results centers paires assignments

centers_dists <- as.matrix( pdist( t(clusters_centres) , t(original_centres) ) )
print(centers_dists)
assign <- solve_LSAP(centers_dists, maximum = FALSE)
print(assign)

########################################
########################################
########################################
########################################

centers_err <- c(1:K)*0

for( d in (1:DIM)){
    unidim_dists <- as.matrix( pdist( t(clusters_centres[,d]) , t(original_centres[,d]) ) )
    print(unidim_dists)
    print(unidim_dists[cbind(seq_along(assign),assign)])
    centers_err[ ] <- centers_err[ ] +  ( unidim_dists[cbind(seq_along(assign),assign)]  / clusters_centres[,d] ) [ ]
}

centers_err <- centers_err / DIM

########################################
########################################
########################################


kmlio_err <- sum(centers_err) / K

dataplus <- data.frame(kmlio_err)

write.table(dataplus,
          paste(result,"log_skip_chunk.csv"),
          quote = FALSE,
          col.names = FALSE,
          row.names = FALSE,
          sep = ",",
          append = T)

# kmlio_err(dists,N,M,K,DIM,SEP,skip_chunks)
