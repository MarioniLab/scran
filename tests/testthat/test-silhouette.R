# Tests the clusterSilhouette() function.
# library(testthat); library(scran); source('setup.R'); source('test-silhouette.R')

set.seed(80000)
test_that('clusterSilhouette yields sensible output for pure clusters', {
    clusters <- rep(seq_len(10), each=100)
    y <- matrix(clusters, nrow=50, ncol=length(clusters), byrow=TRUE)

    out <- clusterSilhouette(y, clusters)
    expect_identical(nrow(out), ncol(y))
    expect_true(all(out$width == 1))
    expect_true(all(clusters != out$other))

    # Throwing in some jitter.
    y <- jitter(y)
    out <- clusterSilhouette(y, clusters)
    expect_true(all(out$width >= 0.5))
    expect_true(all(clusters != out$other))
})

test_that('clusterSilhouette yields correct output for perfectly randomized clusters', {
    clusters <- rep(1:5, each=10)
    y0 <- matrix(rnorm(100), ncol=10)
    y <- cbind(y0, y0, y0, y0, y0)    

    out <- clusterSilhouette(y, clusters)
    expect_identical(nrow(out), ncol(y))
    expect_true(all(out$width == 0))
    expect_true(all(clusters != out$other))
})

set.seed(70002)
test_that('clusterSilhouette handles SCEs correctly', {
    library(scuttle)
    sce <- mockSCE()
    sce <- logNormCounts(sce)
    g <- buildSNNGraph(sce)
    clusters <- igraph::cluster_walktrap(g)$membership

    ref <- clusterSilhouette(logcounts(sce), clusters)
    expect_equal(ref, clusterSilhouette(sce, clusters))

    # Works with base SE's.
    expect_equal(ref, clusterSilhouette(as(sce, "SummarizedExperiment"), clusters))

    # Also deals with reducedDims.
    reducedDim(sce, "stuff") <- matrix(rnorm(10*ncol(sce)), ncol=10)
    expect_equal(
        out <- clusterSilhouette(sce, clusters, use.dimred="stuff"),
        clusterSilhouette(reducedDim(sce, "stuff"), clusters, transposed=TRUE)
    )
    expect_equal(nrow(out), ncol(sce))

    # Works with the supplied clusters.
    colLabels(sce) <- clusters
    expect_equal(ref, clusterSilhouette(sce))
})
