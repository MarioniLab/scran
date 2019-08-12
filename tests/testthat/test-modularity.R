# Tests the clusterModularity function.
# library(scran); library(testthat); source("test-modularity.R")

set.seed(20001)
test_that("clusterModularity computes the correct values", {
    exprs <- matrix(rnorm(100000), ncol=100)
    g <- buildSNNGraph(exprs)

    random <- sample(5, ncol(exprs), replace=TRUE)
    out <- clusterModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random, weight=igraph::E(g)$weight))
    expect_equal(sum(out), 0)

    # Some basic checks on the expected values.
    out <- clusterModularity(g, random, get.weights=TRUE)
    expect_equal(sum(out$observed), sum(out$expected))
    expect_equal(sum(out$observed), sum(igraph::E(g)$weight))
})

set.seed(20001)
test_that("clusterModularity computes correct values in a more realistic scenario", {
    exprs <- matrix(rnorm(100000), ncol=100)
    exprs[1:100,1:20] <- rnorm(2000, 2)
    exprs[101:200,21:50] <-rnorm(3000, 2)
    exprs[201:300,51:90] <-rnorm(4000, 2)
    g <- buildSNNGraph(exprs, k=10)

    actual <- igraph::cluster_walktrap(g)
    out <- clusterModularity(g, actual$membership) 
    expect_equal(sum(diag(out)), igraph::modularity(g, actual$membership, weight=igraph::E(g)$weight))
    expect_equal(sum(out), 0)

    # Some basic checks on the expected values.
    out <- clusterModularity(g, actual$membership, get.weights=TRUE)
    expect_equal(sum(out$observed), sum(out$expected))
    expect_equal(sum(out$observed), sum(igraph::E(g)$weight))
})

set.seed(20002)
test_that("clusterModularity handles unweighted graphs correctly", {
    exprs <- matrix(rnorm(100000), ncol=100)
    g <- buildKNNGraph(exprs)

    random <- sample(5, ncol(exprs), replace=TRUE)
    out <- clusterModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random))
    expect_equal(sum(out), 0)

    # Some basic checks on the expected values.
    out <- clusterModularity(g, random, get.weights=TRUE)
    expect_equal(sum(out$observed), sum(out$expected))
    expect_equal(sum(out$observed), length(igraph::E(g)))
})

set.seed(20003)
test_that("clusterModularity handles directed graphs correctly", {
    exprs <- matrix(rnorm(100000), ncol=100)
    g <- buildKNNGraph(exprs, directed=TRUE)

    random <- sample(5, ncol(exprs), replace=TRUE)
    out <- clusterModularity(g, random)
    ref <- clusterModularity(igraph::as.undirected(g, mode="each"), random) 
    expect_identical(out, ref)
})

set.seed(20003)
test_that("clusterModularity handles self-loops correctly", {
    exprs <- matrix(rnorm(100000), ncol=100)
    g <- buildSNNGraph(exprs)

    g <- igraph::add_edges(g, rep(ncol(exprs), each=2), weight=10)
    random <- sample(5, ncol(exprs), replace=TRUE)
    out <- clusterModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random, weight=igraph::E(g)$weight))

    # Works for unweighted graphs.
    exprs <- matrix(rnorm(100000), ncol=100)
    g <- buildKNNGraph(exprs)

    g <- igraph::add_edges(g, rep(ncol(exprs), each=2))
    random <- sample(5, ncol(exprs), replace=TRUE)
    out <- clusterModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random))
})
