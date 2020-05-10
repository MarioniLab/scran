# Tests the clusterKNNGraph functions and relatives.
# library(testthat); library(scran); source("setup.R"); source("test-cluster-snn.R")

set.seed(1000)
library(scuttle)
sce <- mockSCE(ncells=500, ngenes=1000)
sce <- logNormCounts(sce)

test_that("clusterSNNGraph works correctly for standard applications", {
    clusters <- clusterSNNGraph(sce)
    expect_identical(length(clusters), ncol(sce))

    clusters2 <- clusterSNNGraph(logcounts(sce))
    expect_identical(clusters, clusters2)

    # Works with low-dimensional arguments:
    reducedDim(sce, "stuff") <- matrix(rnorm(10*ncol(sce)), ncol=10)
    clusters3a <- clusterSNNGraph(sce, use.dimred="stuff")
    clusters3b <- clusterSNNGraph(reducedDim(sce, "stuff"), transposed=TRUE)
    expect_identical(clusters3a, clusters3b)

    # Passes arguments down:
    clusters4 <- clusterSNNGraph(sce, k=5)
    expect_false(identical(clusters, clusters4))

    clusters5 <- clusterSNNGraph(sce, clusterFUN=igraph::cluster_louvain)
    expect_false(identical(clusters, clusters5))
})

test_that("clusterSNNGraph works correctly for k-means", {
    clusters <- clusterSNNGraph(sce, use.kmeans=TRUE)
    expect_identical(length(clusters), ncol(sce))

    # Respects the seed:
    set.seed(10)
    clusters2a <- clusterSNNGraph(sce, use.kmeans=TRUE, kmeans.centers=100)
    set.seed(10)
    clusters2b <- clusterSNNGraph(sce, use.kmeans=TRUE, kmeans.centers=100)
    expect_identical(clusters2a, clusters2b)

    # Fetches more details when requested.
    set.seed(10)
    clusters2c <- clusterSNNGraph(sce, use.kmeans=TRUE, kmeans.centers=100, full.stats=TRUE)
    expect_identical(clusters2a, clusters2c$igraph)
    expect_s3_class(metadata(clusters2c)$graph, "igraph")
    expect_identical(clusters2c$igraph, factor(metadata(clusters2c)$membership[clusters2c$kmeans]))
})

test_that("clusterKNNGraph works correctly for standard applications", {
    clusters <- clusterKNNGraph(sce)
    expect_identical(length(clusters), ncol(sce))

    clusters2 <- clusterKNNGraph(logcounts(sce))
    expect_identical(clusters, clusters2)

    # Works with low-dimensional arguments:
    reducedDim(sce, "stuff") <- matrix(rnorm(10*ncol(sce)), ncol=10)
    clusters3a <- clusterKNNGraph(sce, use.dimred="stuff")
    clusters3b <- clusterKNNGraph(reducedDim(sce, "stuff"), transposed=TRUE)
    expect_identical(clusters3a, clusters3b)

    # Passes arguments down:
    clusters4 <- clusterKNNGraph(sce, k=5)
    expect_false(identical(clusters, clusters4))

    clusters5 <- clusterKNNGraph(sce, clusterFUN=igraph::cluster_louvain)
    expect_false(identical(clusters, clusters5))
})

test_that("clusterKNNGraph works correctly for k-means", {
    clusters <- clusterKNNGraph(sce, use.kmeans=TRUE)
    expect_identical(length(clusters), ncol(sce))

    # Respects the seed:
    set.seed(20)
    clusters2a <- clusterKNNGraph(sce, use.kmeans=TRUE, kmeans.centers=100)
    set.seed(20)
    clusters2b <- clusterKNNGraph(sce, use.kmeans=TRUE, kmeans.centers=100)
    expect_identical(clusters2a, clusters2b)

    # Fetches more details when requested.
    set.seed(20)
    clusters2c <- clusterKNNGraph(sce, use.kmeans=TRUE, kmeans.centers=100, full.stats=TRUE)
    expect_identical(clusters2a, clusters2c$igraph)
    expect_s3_class(metadata(clusters2c)$graph, "igraph")
    expect_identical(clusters2c$igraph, factor(metadata(clusters2c)$membership[clusters2c$kmeans]))
})
