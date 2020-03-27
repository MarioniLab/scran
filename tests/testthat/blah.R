source("setup.R")

set.seed(80000)
ncells <- 200
ngenes <- 150
means <- runif(ngenes, 0, 5)
X <- matrix(rpois(ngenes*ncells, lambda=means), ncol=ncells, nrow=ngenes)

library(BiocParallel)
clust <- kmeans(t(X), centers=3)
clusters <- as.factor(clust$cluster)
bplapply(list(X, X), FUN=scater::sumCountsAcrossCells, ids=clusters, BPPARAM=SnowParam(2))

