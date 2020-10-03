# This tests that various functions are applicable with alternative matrix representations.
# library(scran); library(testthat); source("test-altrep.R")

set.seed(99999)
library(Matrix)
X <- as(matrix(rpois(100000, lambda=1), ncol=100), "dgCMatrix")
X_ <- as.matrix(X)

library(HDF5Array)
Y <- as(matrix(rpois(100000, lambda=5), ncol=100), "HDF5Array")
Y_ <- as.matrix(Y)

test_that("cyclone runs properly", {
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    rownames(X) <- rownames(X_) <- sample(mm.pairs$G1[,1], nrow(X))
    rownames(Y) <- rownames(Y_) <- sample(mm.pairs$G1[,1], nrow(Y))

    set.seed(100)
    assignments1 <- cyclone(X[,1:10], mm.pairs)
    set.seed(100)
    assignments2 <- cyclone(X_[,1:10], mm.pairs)
    expect_identical(assignments1, assignments2)

    set.seed(100)
    assignments1 <- cyclone(Y[,1:10], mm.pairs)
    set.seed(100)
    assignments2 <- cyclone(Y_[,1:10], mm.pairs)
    expect_identical(assignments1, assignments2)
})

test_that("Variance estimation runs properly", {
    dec1 <- modelGeneVar(Y)
    dec2 <- modelGeneVar(Y_)
    expect_equal(dec1, dec2)

    dec1 <- modelGeneCV2(Y)
    dec2 <- modelGeneCV2(Y_)
    expect_equal(dec1, dec2)
})

test_that("correlatePairs runs properly", {
    set.seed(1000)
    null <- correlateNull(ncol(X), iters=1e6)

    set.seed(100) 
    ref <- correlatePairs(X_[1:10,], null.dist=null)
    set.seed(100) 
    alt <- correlatePairs(X[1:10,], null.dist=null)
    expect_equal(ref, alt)

    set.seed(200) 
    ref <- correlatePairs(Y_[20:50,], null.dist=null)
    set.seed(200) 
    alt <- correlatePairs(Y[20:50,], null.dist=null)
    expect_equal(ref, alt)
})

test_that("buildSNNGraph with irlba runs properly on sparse matrices", {
    set.seed(1000)
    g1 <- buildSNNGraph(X, BSPARAM=BiocSingular::IrlbaParam())
    set.seed(1000)
    g2 <- buildSNNGraph(X_, BSPARAM=BiocSingular::IrlbaParam())
    expect_identical(g1[], g2[])

    set.seed(100)
    g1 <- buildSNNGraph(X, d=10, BSPARAM=BiocSingular::IrlbaParam())
    set.seed(100)
    g2 <- buildSNNGraph(X_, d=10, BSPARAM=BiocSingular::IrlbaParam())
    expect_identical(g1[], g2[])
})

test_that("findMarkers and overlapExprs work properly", {
    groups <- sample(2, ncol(X), replace=TRUE)
    expect_equal(findMarkers(Y, groups), findMarkers(Y_, groups))
})
