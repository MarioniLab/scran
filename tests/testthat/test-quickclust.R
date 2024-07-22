# This tests out the quickCluster method in scran.
# require(scran); require(testthat); source("setup.R"); source("test-quickclust.R")

set.seed(30000)
ncells <- 700
ngenes <- 1000

dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0
dummy[601:900,known.clusters==3L] <- 0

expect_nested <- function(known, observed) 
# Sometimes the clustering breaks it up into smaller clusters; that's fine,
# as long as those smaller ones are nested within the known clusters.
{
    expect_true(all(lengths(lapply(split(known, observed), unique))==1L))
}

test_that("quickCluster works in the simple case", {
    emp.clusters <- quickCluster(dummy, use.ranks=FALSE)
    expect_nested(known.clusters, emp.clusters)
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

set.seed(30000101)
test_that("use.ranks=TRUE generates the correct ScaledMatrix", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    ref <- scran:::.create_rank_matrix(mat, deferred=FALSE)
    def <- scran:::.create_rank_matrix(mat, deferred=TRUE)
    expect_s4_class(def, "ScaledMatrix")
    expect_equivalent(ref, as.matrix(def))

    # Same results from the two options in quickCluster() itself.
    ref <- quickCluster(mat, use.ranks=TRUE, d=NA, method="hclust")
    bspar <- BiocSingular::ExactParam(deferred=TRUE)
    out <- quickCluster(as(mat, "dgCMatrix"), use.ranks=TRUE, d=min(dim(mat)), method="hclust", BSPARAM=bspar)
    expect_identical(ref, out)

    # Set low 'k' to avoid inconsistencies caused by tied neighbors.
    ref <- quickCluster(mat, use.ranks=TRUE, d=NA, method="igraph", k=2)
    out <- quickCluster(as(mat, "dgCMatrix"), use.ranks=TRUE, d=min(dim(mat)), method="igraph", k=2, BSPARAM=bspar)
    expect_identical(ref, out)
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
    expect_identical(obs, quickCluster(mat, min.size=50, use.ranks=TRUE, 
        subset.row=scuttle::calculateAverage(mat) >= 5))

    # 'min.mean' should have no effect when use.ranks=FALSE.
    obs <- quickCluster(mat, min.size=50, min.mean=5, use.ranks=FALSE)
    expect_identical(obs, quickCluster(mat, min.size=50, min.mean=1, use.ranks=FALSE)) 
})

set.seed(3000012)
test_that("quickCluster functions correctly with blocking", {
    skip_on_os('windows') # 32-bit failure. Who knows, man. Who knows.

    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    block <- sample(3, ncol(mat), replace=TRUE)

    # Using 'hclust' to avoid problems with tied ranks and igraph.
    obs <- quickCluster(mat, min.size=10, block=block, method="hclust", use.ranks=FALSE)

    collected <- numeric(ncol(mat))
    last <- 0L
    for (x in sort(unique(block))) {
        chosen <- block==x
        current <- quickCluster(mat[,chosen], min.size=10, method="hclust", use.ranks=FALSE)
        collected[chosen] <- as.integer(current) + last
        last <- last + nlevels(current)
    }
    expect_identical(obs, factor(collected))

    # Should behave properly with NULL or single-level.
    ref <- quickCluster(mat, min.size=10, block=NULL, method="hclust", use.ranks=FALSE)
    obs <- quickCluster(mat, min.size=10, block=integer(ncol(mat)), method="hclust", use.ranks=FALSE)
    expect_identical(ref, obs)

    # Should avoid problems with multiple BPPARAM specifications.
    ref <- quickCluster(mat, min.size=10, block=block, method="hclust", use.ranks=FALSE)
    obs <- quickCluster(mat, min.size=10, block=block, method="hclust", use.ranks=FALSE, block.BPPARAM=safeBPParam(2))
    expect_identical(obs, ref)
})

set.seed(3000013)
test_that("quickCluster's calls to min.size in dynamic tree cut are respected", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    obs <- scran:::.quick_cluster(mat, min.size=50, method="hclust", d=NA, use.ranks=FALSE)

    ref <- scuttle::normalizeCounts(mat)
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

    forced <- quickCluster(dummy, method="hclust", use.ranks=FALSE, min.size=100) 
    tab <- table(forced)
    nonzero <- setdiff(names(tab), "0")
    expect_true(all(tab[nonzero] >= 70L))
})

set.seed(300002)
test_that("quickCluster with igraph works correctly", {
    k <- 10
    mat <- matrix(rnorm(200000, mean=20), nrow=400)
    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, d=50, use.ranks=FALSE)

    ref <- scuttle::normalizeCounts(mat)
    snn <- buildSNNGraph(ref, k=k, d=50)
    out <- igraph::cluster_walktrap(snn)
    expect_identical(factor(out$membership), obs)

    obs <- quickCluster(mat, min.size=0, method="igraph", k=k, d=50, graph.fun=igraph::cluster_fast_greedy, 
        use.ranks=FALSE)
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
    # These tests are surprisingly fragile for use.rank=TRUE, as findKNN in
    # rank space is liable to find lots of tied distances.  This results in
    # arbitrary choices and ordering of neighbors, which can differ between
    # seeds and machines (depending on precision).  Hence we need to make sure
    # that there are no ties, by supplying enough dimensions with no tied
    # ranks.
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

    # Checking correct subsetting.
    subset.row <- 1:25*2
    set.seed(0) # set.seed() for consistent tie handling by BiocNeighbors.
    clust1 <- quickCluster(X, use.ranks=FALSE, subset.row=subset.row) 
    set.seed(0)
    clust2 <- quickCluster(counts(X)[subset.row,], use.ranks=FALSE)
    expect_identical(clust1, clust2)
})

set.seed(200031)
test_that("quickCluster works on alternative matrices", {
    # Testing with ranks.
    sparse <- abs(Matrix::rsparsematrix(ngenes, ncells, density=0.1))
    out <- quickCluster(sparse, min.mean=0, use.ranks=TRUE)
    ref <- quickCluster(as.matrix(sparse), min.mean=0, use.ranks=TRUE)
    expect_identical(out, ref)

    # Testing without ranks. Note that, if seed is 20003, this results in an
    # extremely unfortunate discrepancy caused by numerical imprecision in the log calculation,
    # such that bio > 0 (barely!) in one run and bio == 0 in another run.
    library(HDF5Array)
    dummy <- as(matrix(rpois(50000, lambda=5), nrow=50), "HDF5Array")
    out <- quickCluster(dummy, min.mean=0, use.ranks=FALSE)
    ref <- quickCluster(as.matrix(dummy), min.mean=0, use.ranks=FALSE)
    expect_identical(out, ref)
})
