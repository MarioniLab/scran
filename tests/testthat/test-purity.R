# Tests the clusterPurity() function.
# library(testthat); library(scran); source('setup.R'); source('test-purity.R')

set.seed(70000)
test_that('clusterPurity yields correct output for pure clusters', {
    clusters <- rep(seq_len(10), each=100)
    y <- matrix(clusters, nrow=50, ncol=length(clusters), byrow=TRUE)
    y <- jitter(y)

    out <- clusterPurity(y, clusters, pseudo.count=0)
    expect_true(all(out==1))
    expect_identical(length(out), ncol(y))
})

set.seed(70001)
test_that('clusterPurity yields correct output for compromised clusters', {
    y <- matrix(seq_len(100), nrow=50, ncol=100, byrow=TRUE)
    y <- jitter(y)
    clusters <- rep(1:2, c(1, 99))
    out <- clusterPurity(y, clusters, pseudo.count=0)

    expect_true(out[1] <= 0.05)
    expect_true(all(out[-1] > 0.9))
    expect_false(is.unsorted(out))
})

set.seed(700011)
test_that('clusterPurity responds to the pseudo-count', {
    y <- matrix(rnorm(10000), nrow=50)
    clusters <- sample(1:5, ncol(y), replace=TRUE)

    out1 <- clusterPurity(y, clusters, pseudo.count=0)
    out2 <- clusterPurity(y, clusters, pseudo.count=1)

    # Greater shrinkage towards 0.5.
    diff1 <- out1 - 0.5
    diff2 <- out2 - 0.5
    expect_true(all(sign(diff1)==sign(diff2)))
    expect_true(all(sign(diff1) * diff1 > sign(diff2) * diff2))
})

set.seed(70002)
test_that('clusterPurity handles SCEs correctly', {
    library(scater)
    sce <- mockSCE()
    sce <- logNormCounts(sce)
    g <- buildSNNGraph(sce)
    clusters <- igraph::cluster_walktrap(g)$membership

    expect_identical(
        clusterPurity(logcounts(sce), clusters),
        clusterPurity(sce, clusters)
    )

    # Handles subsetting correctly.
    expect_identical(
        clusterPurity(sce, clusters, subset.row=1:10),
        clusterPurity(sce[1:10,], clusters)
    )

    # Also deals with reducedDims.
    sce <- runPCA(sce)
    expect_identical(
        out <- clusterPurity(sce, clusters, dimred="PCA"),
        clusterPurity(reducedDim(sce), clusters, transposed=TRUE)
    )
    expect_identical(length(out), ncol(sce))
})
