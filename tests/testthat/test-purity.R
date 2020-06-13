# Tests the clusterPurity() function.
# library(testthat); library(scran); source('setup.R'); source('test-purity.R')

set.seed(70000)
test_that('clusterPurity yields correct output for pure clusters', {
    clusters <- rep(seq_len(10), each=100)
    y <- matrix(clusters, nrow=50, ncol=length(clusters), byrow=TRUE)
    y <- jitter(y)

    out <- clusterPurity(y, clusters)
    expect_identical(nrow(out), ncol(y))

    expect_true(all(out$purity==1))
    expect_true(all(out$maximum==clusters))
})

set.seed(70001)
test_that('clusterPurity yields correct output for compromised clusters', {
    y <- matrix(rep(0:1, each=100), nrow=50, ncol=200, byrow=TRUE)
    y <- jitter(y)

    clusters <- rep(1:2, each=100)
    clusters[1] <- 2

    out <- clusterPurity(y, clusters)

    expect_true(out$purity[1] <= 0.05)
    expect_true(all(out$purity[-1] > 0.9))
    expect_identical(out$maximum[1], 1)
    expect_identical(out$maximum[-1], clusters[-1])
})

set.seed(700001)
test_that('clusterPurity handles subsetting/transposition correctly', {
    y <- matrix(rnorm(10000), nrow=50)
    clusters <- kmeans(t(y), 4)$cluster

    expect_equal(
        clusterPurity(y, clusters, subset.row=1:10),
        clusterPurity(y[1:10,], clusters)
    )

    # Has no effect when transposed.
    expect_equal(
        clusterPurity(t(y), clusters, subset.row=1:10, transposed=TRUE),
        clusterPurity(t(y), clusters, transposed=TRUE)
    )
})


set.seed(700011)
test_that("clusterPurity handles the weighting correctly", {
    # Creating a bulk of points.
    y0 <- matrix(rnorm(10000), nrow=50)
    y1 <- matrix(1000, nrow=50, ncol=10)
    y2 <- matrix(1000, nrow=50, ncol=10)

    y <- cbind(y1, y2, y0)
    clusters <- rep(1:3, c(ncol(y1), ncol(y2), ncol(y0)))
    out1 <- clusterPurity(y, clusters)$purity
    expect_true(all(abs(out1[1:20]-0.5) < 1e-8))
    expect_true(all(out1[-(1:20)]==1))

    # Unaffected by changes in the number of cells. Technically this should be
    # identical(), but the order of cells returned by findNeighbors() changes,
    # and this results in numeric precision changes due to order of addition of
    # double-precision values, resulting in very slightly different output.
    sub <- c(1:5, ncol(y1) + seq_len(ncol(y2) + ncol(y0)))
    out2 <- clusterPurity(y[,sub], clusters[sub])$purity
    expect_equal(out1[sub], out2)

    sub <- c(1, ncol(y1) + seq_len(ncol(y2) + ncol(y0)))
    out2 <- clusterPurity(y[,sub], clusters[sub])$purity
    expect_equal(out1[sub], out2)
})

set.seed(700012)
test_that("clusterPurity handles other weighting options", {
    # Creating a bulk of points.
    y0 <- matrix(rnorm(10000), nrow=50)
    y1 <- matrix(1000, nrow=50, ncol=10)
    y2 <- matrix(1000, nrow=50, ncol=10)

    y <- cbind(y1, y2, y0)
    clusters <- rep(1:3, c(ncol(y1), ncol(y2), ncol(y0)))

    # Turning off weighting has no effect for balanced clusters.
    out1 <- clusterPurity(y, clusters, weighted=FALSE)$purity
    expect_true(all(abs(out1[1:20]-0.5) < 1e-8))
    expect_true(all(out1[-(1:20)]==1))

    # Turning off weighting has some effect for non-balanced clusters.
    sub <- c(1:5, ncol(y1) + seq_len(ncol(y2) + ncol(y0)))
    out2 <- clusterPurity(y[,sub], clusters[sub], weighted=FALSE)$purity
    expect_true(all(abs(out2[1:5]-1/3) < 1e-8))
    expect_true(all(abs(out2[6:15]-2/3) < 1e-8))
    expect_true(all(out2[-(1:20)]==1))

    # We can replace it with our own weighting to restore the balance.
    out3 <- clusterPurity(y[,sub], clusters[sub], weighted=rep(c(2, 1, 1), c(5, ncol(y2), ncol(y0))))$purity
    ref <- clusterPurity(y[,sub], clusters[sub])$purity
    expect_equal(out3, ref)
})

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
