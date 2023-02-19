
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) # alternatively, this also loads %>%
library(mclust)

folder <- "generator/CM5,8M_1000MO_SEP0,3/"

# ari_classes <- read.table(paste(folder,"datapoints_1.mem",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)

# ari_clusters <- read.table("./results/clusters.csv", header = FALSE,  sep = "\t", stringsAsFactors = FALSE)

# mc_ari <- adjustedRandIndex(ari_clusters$V1, ari_classes$V1)
# print(mc_ari)
mc_ari <- -1
# ari_ct <- table(ari_classes$V1, ari_clusters$V1)
# print(ari_ct)

#############################################
#             clusters centers
#############################################

init_method <- 5
N <- 500
M <- 100
K <- 10
V <- 10
SEP <- 0.3

original_centres <- read.table(paste(folder,"original_centres_norm.csv",sep=""),header = FALSE,sep = "\t",stringsAsFactors = FALSE)

clusters_centres <- read.table("results/centers.csv", header = FALSE,  sep = "\t", stringsAsFactors = FALSE)

library(pdist)
dists <- pdist(t(clusters_centres[,-11]), t(original_centres))
as.matrix(dists)

library("sos")
#findFn("colMins")
library("matrixStats")
min_dists <- colMins(as.matrix(dists))

dataplus <- data.frame(init_method, N, M, N/M, K, V, SEP, mc_ari,sum(min_dists)/K,-1,-1,-2,t(min_dists))

write.table(dataplus,
          "generator/clusters_eval_log.csv",
          quote = FALSE,
          col.names = FALSE,
          row.names = FALSE,
          sep = "\t",
          append = T) # nolint

#ari_cc <- ARI(ari_clusters$V1, ari_classes$V1)
#all.equal(mc_ari, unclass(ari_cc), check.attributes = FALSE)
