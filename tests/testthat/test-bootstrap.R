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
    expect_true(all(output[upper.tri(output, diag=FALSE)] <= 0.1))
    expect_true(all(diag(output) > 0.5))

    # Continues to work if vector is a character or factor.
    output <- bootstrapCluster(dummy, FUN=function(x) { c("X", "Y", "Z")[kmeans(t(log10(x+1)), 3)$cluster] })
    expect_true(all(output[upper.tri(output, diag=FALSE)] <= 0.1))
    expect_true(all(diag(output) > 0.5))

    output <- bootstrapCluster(dummy, FUN=function(x) { factor(kmeans(t(log10(x+1)), 3)$cluster) })
    expect_true(all(output[upper.tri(output, diag=FALSE)] <= 0.1))
    expect_true(all(diag(output) > 0.5))
})

set.seed(500002)
test_that("bootstrapCluster works correctly with poor separation", {
    dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
    output <- bootstrapCluster(dummy, FUN=function(x) { kmeans(t(log10(x+1)), 3)$cluster })
    expect_true(all(output[upper.tri(output, diag=TRUE)] > 0.2))
    expect_true(all(diag(output) < 0.8))
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
