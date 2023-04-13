
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) # alternatively, this also loads %>%
library(mclust)
library(clue)

library(pdist)
library("sos")
library("matrixStats")


kmlio_err <- function(centers_dists,N,M,K,DIM,SEP,skip_chunks) {
    # print(centers_dists)
    assign <- solve_LSAP(centers_dists, maximum = FALSE)
    print(assign)
    ## To get the optimal value (for now):
    dist_err <- sum(centers_dists[cbind(seq_along(assign),assign)]) / K

    dataplus <- data.frame(N, M, N/M, K, DIM, SEP,skip_chunks,dist_err)

    write.table(dataplus,
          "generator/kmlio_cluster_evaluation.csv",
          quote = FALSE,
          col.names = FALSE,
          row.names = FALSE,
          sep = "\t",
          append = T)

    # done
}

#############################################
#             clusters centers
#############################################

folder <- "generator/10D/13421800N/SEP0.3/"

N <- 13421800
M <- 3355450    #6710900
K <- 10
DIM <- 10
V <- 10
SEP <- 0.3
skip_chunks <-0

original_centres <- read.table(paste(folder,"real_centers_norm.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)
clusters_centres <- read.table("results/result_centers.csv", header = FALSE,  sep = "\t", stringsAsFactors = FALSE)

dists <- as.matrix(pdist(t(clusters_centres[,-11]), t(original_centres)))
print(dists)

kmlio_err(dists,N,M,K,DIM,SEP,skip_chunks)
