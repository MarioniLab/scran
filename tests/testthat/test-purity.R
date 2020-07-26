# Tests the clusterPurity() function.
# library(testthat); library(scran); source('setup.R'); source('test-purity.R')

set.seed(70002)
test_that('clusterPurity handles SCEs correctly', {
    library(scuttle)
    sce <- mockSCE()
    sce <- logNormCounts(sce)
    g <- buildSNNGraph(sce)
    clusters <- igraph::cluster_walktrap(g)$membership

    ref <- clusterPurity(logcounts(sce), clusters)
    expect_equal(ref, clusterPurity(sce, clusters))

    # Works with base SE's.
    expect_equal(ref, clusterPurity(as(sce, "SummarizedExperiment"), clusters))

    # Also deals with reducedDims.
    reducedDim(sce, "stuff") <- matrix(rnorm(10*ncol(sce)), ncol=10)
    expect_equal(
        out <- clusterPurity(sce, clusters, use.dimred="stuff"),
        clusterPurity(reducedDim(sce, "stuff"), clusters, transposed=TRUE)
    )
    expect_equal(nrow(out), ncol(sce))

    # Works with the supplied clusters.
    colLabels(sce) <- clusters
    expect_equal(ref, clusterPurity(sce))
})
