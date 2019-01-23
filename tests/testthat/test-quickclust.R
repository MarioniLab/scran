# This tests out the quickCluster method in scran.
# require(scran); require(testthat); source("test-quickclust.R")

set.seed(30000)
ncells <- 700
ngenes <- 1000

dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0
dummy[601:900,known.clusters==3L] <- 0

test_that("quickCluster works in the simple case", {
    emp.clusters <- quickCluster(dummy, use.ranks=FALSE)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
    emp.clusters <- quickCluster(dummy, d=20, use.ranks=FALSE)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)

    # Behaves after rows are shuffled.
    shuffled <- sample(nrow(dummy))
    shuf.clusters <- quickCluster(dummy, subset.row=shuffled, use.ranks=FALSE)
    expect_true(length(unique(paste0(shuf.clusters, emp.clusters)))==3L)

    # Behaves with use.ranks=TRUE.
    emp.clusters <- quickCluster(dummy, use.ranks=TRUE)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
    emp.clusters <- quickCluster(dummy, use.ranks=TRUE, d=20)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)

    # Behaves when PCA is turned off.
    emp.clusters <- quickCluster(dummy, d=NA, use.ranks=FALSE)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
    emp.clusters <- quickCluster(dummy, d=NA, use.ranks=TRUE)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
})

set.seed(300001)
test_that("quickCluster with use.ranks=TRUE is consistent with clustering on correlations", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    obs <- quickCluster(mat, use.ranks=TRUE, method="hclust", d=NA, min.size=20)

    refM <- sqrt(0.5*(1 - cor(mat, method="spearman")))
    distM <- as.dist(refM)
    obsM <- dist(scaledColRanks(mat, transposed=TRUE))
    expect_equal(as.matrix(obsM), as.matrix(distM))

    htree <- hclust(distM, method='ward.D2')
    clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=20, distM=refM, verbose=0))
    expect_identical(clusters, as.integer(obs)) # this can complain if unassigned, as 0 becomes 1 in as.integer().
})

set.seed(3000011)
test_that("quickCluster functions correctly with subsetting", {
    mat <- matrix(rpois(20000, lambda=1:100), nrow=100)

    # Works properly with subsetting.
    subset.row <- sample(nrow(mat), nrow(mat)/2)
    obs <- quickCluster(mat, min.size=50, subset.row=subset.row, use.ranks=FALSE)
    expect_identical(obs, quickCluster(mat[subset.row,], min.size=50, use.ranks=FALSE)) # Checking that it behaves properly.
    expect_false(identical(quickCluster(mat, min.size=50, use.ranks=FALSE), obs)) # It should return different results.

    # Same for ranks (use hclust to avoid issues with tied neighbors in rank space).
    obs <- quickCluster(mat, min.size=50, subset.row=subset.row, use.ranks=TRUE, method="hclust")
    expect_identical(obs, quickCluster(mat[subset.row,], min.size=50, use.ranks=TRUE, method="hclust"))
    expect_false(identical(quickCluster(mat, min.size=50, use.ranks=TRUE, method="hclust"), obs))

    # Handles the mean, but only when use.ranks=TRUE.
    obs <- quickCluster(mat, min.size=50, min.mean=5, use.ranks=TRUE)
    expect_identical(obs, quickCluster(mat, min.size=50, use.ranks=TRUE, subset.row=scater::calcAverage(mat) >= 5))

    obs <- quickCluster(mat, min.size=50, min.mean=5, use.ranks=FALSE)
    expect_identical(obs, quickCluster(mat, min.size=50, min.mean=1, use.ranks=FALSE)) # should not respond.
})

set.seed(3000012)
test_that("quickCluster functions correctly with blocking", {
    # Comparison to a slow manual method
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    block <- sample(3, ncol(mat), replace=TRUE)
    obs <- quickCluster(mat, min.size=10, block=block, use.ranks=FALSE)

    collected <- numeric(ncol(mat))
    last <- 0L
    for (x in sort(unique(block))) {
        chosen <- block==x
        current <- quickCluster(mat[,chosen], min.size=10, use.ranks=FALSE)
        collected[chosen] <- as.integer(current) + last
        last <- last + nlevels(current)
    }
    expect_identical(obs, factor(collected))

    # Should behave properly with NULL or single-level.
    ref <- quickCluster(mat, min.size=10, block=NULL, use.ranks=FALSE)
    obs <- quickCluster(mat, min.size=10, block=integer(ncol(mat)), use.ranks=FALSE)
    expect_identical(ref, obs)

    # Should avoid problems with multiple BPPARAM specifications.
    ref <- quickCluster(mat, min.size=10, block=block, method="igraph", use.ranks=FALSE)
    obs <- quickCluster(mat, min.size=10, block=block, method="igraph", use.ranks=FALSE, block.BPPARAM=BiocParallel::MulticoreParam())
    expect_identical(obs, ref)
})

set.seed(3000013)
test_that("quickCluster's calls to min.size in dynamic tree cut are respected", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    obs <- scran:::.quick_cluster(mat, min.size=50, method="hclust", d=NA, use.ranks=FALSE)

    ref <- scater::normalizeCounts(mat, scater::librarySizeFactors(mat), return_log=TRUE)
    refM <- dist(t(ref))
    htree <- hclust(refM, method="ward.D2")
    clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=as.matrix(refM), verbose=0))
    expect_identical(clusters, as.integer(obs))

    # Forcing min.size to be larger than a population size, checking for a warning upon unassigned cells.
    ncells <- 200
    dummy <- matrix(rpois(ncells*200, lambda=5), nrow=200)
    known.clusters <- sample(3, ncells, replace=TRUE)
    dummy[1:40,known.clusters==1L] <- 0
    dummy[41:80,known.clusters==2L] <- 0
    dummy[81:120,known.clusters==3L] <- 0

    out <- quickCluster(dummy, min.size=0, method="hclust", use.ranks=FALSE)
    expect_identical(length(unique(paste(out, known.clusters))), 3L)

    leftovers <- min(80, sum(known.clusters==3))
    keep <- c(which(known.clusters==1), which(known.clusters==2), which(known.clusters==3)[seq_len(leftovers)]) # force cluster 3 to be unassigned.
    expect_warning(forced <- quickCluster(dummy[,keep], method="hclust", use.ranks=FALSE, min.size=100), 
        sprintf("%i cells were not assigned to any cluster", leftovers))
    expect_identical(as.character(tail(forced, leftovers)), rep("0", leftovers))
})

set.seed(300002)
test_that("quickCluster with igraph works correctly", {
    k <- 10
    mat <- matrix(rnorm(200000, mean=20), nrow=400)
    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, d=50, use.ranks=FALSE)

    ref <- scater::normalizeCounts(mat, scater::librarySizeFactors(mat), return_log=TRUE)
    snn <- buildSNNGraph(ref, k=k, d=50)
    out <- igraph::cluster_walktrap(snn)
    expect_identical(factor(out$membership), obs)

    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, d=50, graph.fun=igraph::cluster_fast_greedy, use.ranks=FALSE)
    out <- igraph::cluster_fast_greedy(snn)
    expect_identical(factor(out$membership), obs)

    # Checking that 'd' is respected, along with other arguments that are passed along.
    obs <- quickCluster(mat, min.size=0, method="igraph", d=20, k=15, use.ranks=FALSE)
    snn <- buildSNNGraph(ref, d=20, k=15)
    out <- igraph::cluster_walktrap(snn)
    expect_identical(factor(out$membership), obs)

    obs <- quickCluster(mat, min.size=0, method="igraph", d=NA, k=5, use.ranks=FALSE)
    snn <- buildSNNGraph(ref, d=NA, k=5)
    out <- igraph::cluster_walktrap(snn)
    expect_identical(factor(out$membership), obs)
})

set.seed(3000021)
test_that("quickCluster with igraph merging works correctly", {
    k <- 10
    mat <- matrix(rnorm(200000, mean=20), nrow=400)
    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, use.ranks=FALSE)

    min.size <- 100
    expect_false(all(table(obs) >= min.size))
    obs2 <- quickCluster(mat, min.size=min.size, method="igraph", k=k, use.ranks=FALSE)
    expect_true(all(table(obs2) >= min.size))

    combined <- paste0(obs, ".", obs2)
    expect_identical(length(unique(combined)), length(unique(obs))) # Confirm that they are nested.
})

set.seed(3000022)
test_that("quickCluster with igraph on ranks works correctly", {
    # NOTE 1: these tests are surprisingly fragile for use.rank=TRUE, as findKNN in rank space is liable to find lots of tied distances.
    # This results in arbitrary choices and ordering of neighbors, which can differ between seeds and machines (depending on precision).
    # Hence we need to make sure that there are no ties, by supplying enough dimensions with no tied ranks.
    #
    # NOTE 2: As a result of the above note, we also need to turn off approximate PCs, to avoid issues with irlba variability.
    # Otherwise we would have to set the seed everytime, and this would mask any issues with NN detection.
    k <- 10
    mat <- matrix(rnorm(200000, mean=20), nrow=400)
    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, use.ranks=FALSE)

    # Checking that there are no ties within the 'k+1'th nearest neighbors for each cell.
    ref <- scaledColRanks(mat)
    all.dist <- as.matrix(dist(t(ref)))
    diag(all.dist) <- Inf # ignore self.
    out <- apply(all.dist, 1, FUN=function(d) { min(diff(sort(d)[seq_len(k+1)])) })
    expect_true(min(out) > 1e-8)

    # Testing igraph mode.
    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, use.ranks=TRUE)
    expect_identical(length(obs), ncol(mat))

    snn <- buildSNNGraph(ref, k=k, d=50)
    out <- igraph::cluster_walktrap(snn)
    expect_identical(factor(out$membership), obs)

    # Passes on other parameters.
    obs <- quickCluster(mat, min.size=0, method="igraph", k=10, d=20, use.ranks=TRUE)
    snn <- buildSNNGraph(ref, k=10, d=20)
    out <- igraph::cluster_walktrap(snn)
    expect_identical(factor(out$membership), obs)
})

test_that("quickCluster fails on silly inputs", {
    dummy <- matrix(rpois(10000, lambda=5), nrow=20)
    expect_error(quickCluster(dummy[0,], use.ranks=FALSE), "need at least 2 points")
    expect_error(quickCluster(dummy[,0], use.ranks=FALSE), "no residual d.f.")

    expect_error(quickCluster(dummy[0,], use.ranks=TRUE), "rank variances of zero detected for a cell")
    expect_error(quickCluster(dummy[,0], use.ranks=TRUE), "a dimension is zero")
    expect_error(quickCluster(dummy[,0], d=NA, use.ranks=TRUE), "fewer cells than the minimum cluster size")
})

set.seed(20002)
test_that("quickCluster works on SingleCellExperiment objects", {
    dummy <- matrix(rpois(50000, lambda=5), nrow=50)
    rownames(dummy) <- paste0("X", seq_len(nrow(dummy)))
    X <- SingleCellExperiment(list(counts=dummy))
    emp.clusters <- quickCluster(X, use.ranks=FALSE)
    expect_identical(emp.clusters, quickCluster(counts(X), use.ranks=FALSE))

    # Checking correct interplay between spike-ins and subset.row.
    isSpike(X, "ERCC") <- 1:20
    expect_identical(quickCluster(X, use.ranks=FALSE), quickCluster(counts(X)[-(1:20),], use.ranks=FALSE))

    subset.row <- 1:25*2
    set.seed(0); clust1 <- quickCluster(X, use.ranks=FALSE, subset.row=subset.row) # set.seed() for consistent tie handling by BiocNeighbors.
    set.seed(0); clust2 <- quickCluster(counts(X)[setdiff(subset.row, 1:20),], use.ranks=FALSE)
    expect_identical(clust1, clust2)

    expect_identical(quickCluster(X, subset.row=subset.row, get.spikes=TRUE, use.ranks=FALSE),
                     quickCluster(counts(X)[subset.row,], use.ranks=FALSE))
})

set.seed(20003)
test_that("quickCluster works on alternative matrices", {
    # Testing with ranks.
    sparse <- abs(Matrix::rsparsematrix(ngenes, ncells, density=0.1))
    out <- quickCluster(sparse, min.mean=0, use.ranks=TRUE)
    ref <- quickCluster(as.matrix(sparse), min.mean=0, use.ranks=TRUE)
    expect_identical(out, ref)

    # Testing without ranks.
    library(HDF5Array)
    dummy <- as(matrix(rpois(50000, lambda=5), nrow=50), "HDF5Array")
    out <- quickCluster(dummy, min.mean=0, use.ranks=FALSE)
    ref <- quickCluster(as.matrix(dummy), min.mean=0, use.ranks=FALSE)
    expect_identical(out, ref)
})
