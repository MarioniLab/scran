# This tests the bootstrapCluster functionality.
# library(testthat); library(scran); source("setup.R"); source("test-bootstrap.R")

set.seed(50000)
ncells <- 700
ngenes <- 1000

set.seed(500001)
test_that("bootstrapCluster works correctly with clear separation", {
    dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
    known.clusters <- sample(3, ncells, replace=TRUE)
    dummy[1:300,known.clusters==1L] <- 0
    dummy[301:600,known.clusters==2L] <- 0
    dummy[601:900,known.clusters==3L] <- 0

    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(t(log10(x+1)), 3)$cluster })
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))

    # Works with the mean.
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(t(log10(x+1)), 3)$cluster }, average="mean")
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))

    # Continues to work if vector is a character or factor.
    output <- bootstrapCluster(dummy, FUN=function(x) { c("X", "Y", "Z")[kmeans(t(log10(x+1)), 3)$cluster] })
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))

    output <- bootstrapCluster(dummy, FUN=function(x) { factor(kmeans(t(log10(x+1)), 3)$cluster) })
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))
})

set.seed(500002)
test_that("bootstrapCluster works correctly with poor separation", {
    dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(t(log10(x+1)), 3)$cluster })

    expect_true(all(output[upper.tri(output, diag=TRUE)] < 0.1))
    expect_true(all(output[upper.tri(output, diag=TRUE)] > -0.1))
    expect_true(all(diag(output) < 0.1))
    expect_true(all(diag(output) > -0.1))

    # Works with the mean.
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(t(log10(x+1)), 3)$cluster }, average="mean")
    expect_true(all(output[upper.tri(output, diag=TRUE)] < 0.1))
    expect_true(all(output[upper.tri(output, diag=TRUE)] > -0.1))
    expect_true(all(diag(output) < 0.1))
    expect_true(all(diag(output) > -0.1))
})

set.seed(500003)
test_that("bootstrapCluster works correctly with transposed inputs.", {
    dummy <- matrix(rnorm(ncells*20), nrow=ncells, ncol=20)
    
    set.seed(10)
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(x, 3)$cluster }, transposed=TRUE)
    set.seed(10)
    ref <- bootstrapCluster(t(dummy), FUN=function(x) { kmeans(t(x), 3)$cluster })

    expect_identical(output, ref)
})

set.seed(500004)
test_that("bootstrapCluster works when some clusters are not in the bootstrap.", {
    dummy <- matrix(rnorm(10), ncol=10)

    # Guaranteed to get missing clusters from resampling.
    output <- bootstrapCluster(dummy, FUN=function(x) { seq_len(ncol(x)) })

    expect_identical(rownames(output), as.character(seq_len(ncol(dummy))))
    expect_true(all(output[upper.tri(output, diag=FALSE)]==0))
})

set.seed(500003)
test_that("bootstrapCluster works with alternative comparison functions", {
    dummy <- matrix(rnorm(ncells*20), nrow=ncells, ncol=20)

    set.seed(10)
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(x, 3)$cluster }, compare=clusterRand)
    expect_true(all(is.na(output[lower.tri(output)])))
    expect_true(all(!is.na(output[!lower.tri(output)])))
})

set.seed(500004)
test_that("other miscellaneous tests for bootstrapCluster", {
    dummy <- matrix(rnorm(ncells*20), nrow=ncells, ncol=20)

    set.seed(20)
    ref <- bootstrapCluster(dummy, FUN=function(x) { kmeans(x, 3)$cluster }, transposed=TRUE)
    set.seed(20)
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(x, 3)$cluster }, transposed=TRUE)
    
    expect_identical(ref, output)
    expect_error(bootstrapCluster(dummy, FUN=function(x) { seq_len(ncol(x)) }, iterations=0), "positive")
})
