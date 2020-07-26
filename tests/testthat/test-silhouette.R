# Tests the clusterSilhouette() function.
# library(testthat); library(scran); source('setup.R'); source('test-silhouette.R')

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
