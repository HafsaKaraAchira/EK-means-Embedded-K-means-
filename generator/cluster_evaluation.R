if(!require("clue")) { 
    install.packages("clue")
    library("clue")
}

if(!require("pdist")) { 
    install.packages("pdist")
    library("pdist")
}

if(!require("sos")) { 
    install.packages("sos")
    library("sos")
}

if(!require("matrixStats")) { 
    install.packages("matrixStats")
    library("matrixStats")
}

library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) # alternatively, this also loads %>%
library(mclust)

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
BETA <- args[12]
EXP_ID <- args[13]

report <- sprintf("SKIP_CHUNKS/reports/SEP%s/%sN_%sM_%sD_%sK_%sL_%sbeta_EXP%s/",SEP,N,M,DIM,K,DMAX,BETA,EXP_ID)
folder <- sprintf("generator/%sD/%sN/SEP%s/",DIM,N,SEP) 


ari_classes <- read.table(paste(folder,"real_classes.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE,blank.lines.skip = TRUE,nrows = as.numeric(N))

ari_clusters <- read.table(paste(report,"result_clusters.csv",sep="") , header = FALSE,  sep = "\t", stringsAsFactors = FALSE,blank.lines.skip = TRUE,nrows = as.numeric(N))

kmlio_ari <- adjustedRandIndex(ari_clusters$V1, ari_classes$V1)
print(kmlio_ari)

ari_ct <- table(ari_classes$V1, ari_clusters$V1)
print(ari_ct)

dataplus <- data.frame(ari_ct)

write.table(
        dataplus,
        paste(report,"ARI_CT.csv",sep=""),
        quote = FALSE,
        col.names = FALSE,
        row.names = FALSE,
        sep = ",",
        append = F
        )

########################################
########################################

original_centres <- read.table(paste(folder,"real_centers_norm.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE,blank.lines.skip = TRUE,nrows = as.numeric(K))[,1:DIM]
clusters_centres <- read.table(paste(report,"result_centers.csv",sep=""), header = FALSE,  sep = "\t", stringsAsFactors = FALSE,blank.lines.skip = TRUE,nrows = as.numeric(K))[,1:DIM]

glimpse(original_centres)
glimpse(clusters_centres)

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
print(centers_err)

for( d in (1:DIM)){
    unidim_dists <- as.matrix( pdist( data.frame(clusters_centres[,d]) , data.frame(original_centres[,d]) ) )
    # print(unidim_dists)
    print("unidimensional distance :\n")
    print(unidim_dists[cbind(seq_along(assign),assign)])
    centers_err[ ] <- centers_err[ ] +  ( unidim_dists[cbind(seq_along(assign),assign)]  / clusters_centres[,d] ) [ ]
    print("centers error sum :\n")
    print(centers_err)
}

centers_err <- centers_err / as.numeric(DIM)
print(centers_err)

########################################
########################################
########################################


kmlio_err <- sum(centers_err) / as.numeric(K)
print(kmlio_err)

dataplus <- data.frame(kmlio_err,kmlio_ari)

write.table(
        dataplus,
        "SKIP_CHUNKS/reports/log_skip_chunk.csv",
        quote = FALSE,
        col.names = FALSE,
        row.names = FALSE,
        sep = ",",
        append = T
        )

# kmlio_err(dists,N,M,K,DIM,SEP,skip_chunks)
