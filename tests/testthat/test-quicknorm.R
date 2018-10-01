# This tests the quickSumFactors function.
# library(testthat); library(scran); source("test-quicknorm.R")

COREFUN <- function(x, indices, distances, ref.cell, min.mean=0, ndist=3) {
    lib.sizes <- colSums(x)
    local.ave <- t(scran:::.compute_tricube_average(t(x), indices, distances, ndist=ndist))
    local.lib <- colSums(local.ave)

    ref.ave <- local.ave[,ref.cell]
    ref.lib <- local.lib[ref.cell]
    ratios <- local.ave/ref.ave

    sf <- numeric(ncol(x))
    for (j in seq_len(ncol(local.ave))) {
        ab <- (local.ave[,j] / local.lib[j] + ref.ave / ref.lib) / 2 * (local.lib[j] + ref.lib)/2
        keep <- ab >= min.mean
        sf[j] <- median(ratios[keep,j], na.rm=TRUE) * lib.sizes[j] / local.lib[j]
    }

    return(list(sf=sf, ref=ref.ave))
}

test_that("C++ code underlying quickSumFactors() works correctly", {
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) # using rgamma() to avoid ties.
    closest <- BiocNeighbors::findKNN(t(x), k=50)

    # No filtering.
    ref <- COREFUN(x, closest$index, closest$distance, 1, 0)
    test <- scran:::.quick_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=NULL)
    expect_equal(ref, test)

    x2 <- x 
    x2[1:10,] <- 0 # handles all-zeroes properly.
    ref <- COREFUN(x2, closest$index, closest$distance, 1, 0)
    test <- scran:::.quick_sum_cpp_wrapper(x2, closest$index, closest$distance, 1L, min.mean=NULL)
    expect_equal(ref, test)

    x2 <- x 
    x2[sample(length(x2), 10000)] <- 0 # handles random zeroes properly.
    ref <- COREFUN(x2, closest$index, closest$distance, 1, 0)
    test <- scran:::.quick_sum_cpp_wrapper(x2, closest$index, closest$distance, 1L, min.mean=NULL)
    expect_equal(ref, test)

    # Plus filtering.
    ref <- COREFUN(x, closest$index, closest$distance, 1, min.mean=1)
    test <- scran:::.quick_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=1)
    expect_equal(ref, test)

    ref <- COREFUN(x, closest$index, closest$distance, 1, min.mean=2)
    test <- scran:::.quick_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=2)
    expect_equal(ref, test)

    # Different reference cell.
    ref <- COREFUN(x, closest$index, closest$distance, 10L, min.mean=0)
    test <- scran:::.quick_sum_cpp_wrapper(x, closest$index, closest$distance, 10L, min.mean=NULL)
    expect_equal(ref, test)

    # Different ndist scaling.
    ref <- COREFUN(x, closest$index, closest$distance, 1L, min.mean=0, ndist=2)
    test <- scran:::.quick_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=NULL, ndist=2)
    expect_equal(ref, test)
})

test_that("quickSumFactors() works correctly in vanilla applications", {
    # An overall run to check that it behaves.
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) 
    out <- quickSumFactors(x)
    expect_equal(length(out), ncells)
    expect_equal(mean(out), 1)

    # A silly run to check that it returns the true size factors exactly.
    true.sf <- runif(ncells)
    x <- outer(rpois(ngenes, lambda=20), true.sf)
    expect_warning(out <- quickSumFactors(x), "tied")
    expect_equal(out, true.sf/mean(true.sf))

    # Another run to check that subpopulations are properly identified.
    ncells <- 600
    ngenes <- 200
    count.sizes <- rnbinom(ncells, mu=100, size=5)
    multiplier <- seq_len(ngenes)/100
    dummy <- outer(multiplier, count.sizes)

    known.clusters <- sample(3, ncells, replace=TRUE) # Most genes (120 out of 200) are DE in at least one cluster.
    dummy[1:40,known.clusters==1L] <- 0
    dummy[41:80,known.clusters==2L] <- 0
    dummy[81:120,known.clusters==3L] <- 0

    expect_warning(out <- quickSumFactors(dummy), "tied")
    expect_equal(out, count.sizes/mean(count.sizes))
})

