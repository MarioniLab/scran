# This tests the simpleSumFactors function.
# library(testthat); library(scran); source("test-simple-norm.R")

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

set.seed(120000)
test_that("C++ code underlying simpleSumFactors() works correctly", {
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) # using rgamma() to avoid ties.
    closest <- BiocNeighbors::findKNN(t(x), k=50)

    # No filtering.
    ref <- COREFUN(x, closest$index, closest$distance, 1, 0)
    test <- scran:::.simple_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=NULL)
    expect_equal(ref, test)

    x2 <- x 
    x2[1:10,] <- 0 # handles all-zeroes properly.
    ref <- COREFUN(x2, closest$index, closest$distance, 1, 0)
    test <- scran:::.simple_sum_cpp_wrapper(x2, closest$index, closest$distance, 1L, min.mean=NULL)
    expect_equal(ref, test)

    x2 <- x 
    x2[sample(length(x2), 10000)] <- 0 # handles random zeroes properly.
    ref <- COREFUN(x2, closest$index, closest$distance, 1, 0)
    test <- scran:::.simple_sum_cpp_wrapper(x2, closest$index, closest$distance, 1L, min.mean=NULL)
    expect_equal(ref, test)

    # Plus filtering.
    ref <- COREFUN(x, closest$index, closest$distance, 1, min.mean=1)
    test <- scran:::.simple_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=1)
    expect_equal(ref, test)

    # Different reference cell.
    ref <- COREFUN(x, closest$index, closest$distance, 10L, min.mean=0)
    test <- scran:::.simple_sum_cpp_wrapper(x, closest$index, closest$distance, 10L, min.mean=NULL)
    expect_equal(ref, test)

    # Different ndist scaling.
    ref <- COREFUN(x, closest$index, closest$distance, 1L, min.mean=0, ndist=2)
    test <- scran:::.simple_sum_cpp_wrapper(x, closest$index, closest$distance, 1L, min.mean=NULL, ndist=2)
    expect_equal(ref, test)
})

set.seed(120001)
test_that("simpleSumFactors() works correctly in vanilla applications", {
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) 
    out <- simpleSumFactors(x)
    expect_equal(length(out), ncells)
    expect_equal(mean(out), 1)
    expect_true(cor(out, colSums(x)) > 0.9) # basically library size normalization.
})

set.seed(120002)
test_that("simpleSumFactors() responds to all parameters", {
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) 

    # First, a positive control, to check that our tolerances are ok.
    ref <- simpleSumFactors(x)
    expect_true(isTRUE(all.equal(ref, simpleSumFactors(x)))) 

    # Check that parameter values are passed down properly and modify the result.
    expect_false(isTRUE(all.equal(ref, simpleSumFactors(x, k=19))))
    expect_false(isTRUE(all.equal(ref, simpleSumFactors(x, min.mean=0.95))))
    expect_false(isTRUE(all.equal(ref, simpleSumFactors(x, BNPARAM=BiocNeighbors::AnnoyParam(ntrees=10)))))

    # Checking that subsetting works as expected.
    expect_equal(simpleSumFactors(x[1:500,]), simpleSumFactors(x, subset.row=1:500))
})

set.seed(120003)
test_that("simpleSumFactors() works correctly with blocks", {
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) 
    out <- simpleSumFactors(x, block=gl(2, 100))
    expect_equal(length(out), ncells)
    expect_equal(mean(out), 1)

    out1 <- simpleSumFactors(x[,1:100])
    out2 <- simpleSumFactors(x[,1:100+100])
    sub1 <- out[1:100]
    expect_equal(out1, sub1/mean(sub1))
    sub2 <- out[1:100+100]
    expect_equal(out2, sub2/mean(sub2))

    # Checking that it behaves when block indices are unordered.
    scramble <- sample(ncells)
    out2 <- simpleSumFactors(x[,scramble], block=gl(2, 100)[scramble])
    expect_equal(out[scramble], out2)

    # Checking that it behaves in parallel.
    out3 <- simpleSumFactors(x, block=gl(2, 100), BPPARAM=MulticoreParam(2))
    expect_equal(out3, out)

    # Checking that rescaling is done correctly.
    x2 <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) 
    x3 <- matrix(rgamma(ngenes*ncells, 2, 2), nrow=ngenes, ncol=ncells) 

    block <- gl(3, ncells)
    out <- simpleSumFactors(cbind(x, x2*2, x3*3), block=block)
    expect_equal(mean(out[block==2])/mean(out[block==1]), 2, tol=0.1)
    expect_equal(mean(out[block==3])/mean(out[block==1]), 3, tol=0.1)
})

set.seed(120004)
test_that("simpleSumFactors() deals with alternative representations", {
    ncells <- 200
    ngenes <- 1000

    library(Matrix)
    y <- rsparsematrix(ngenes, ncells, density=0.1)
    y <- abs(y)
    expect_s4_class(y, "dgCMatrix")
    expect_equal(simpleSumFactors(y, min.mean=0), simpleSumFactors(as.matrix(y), min.mean=0))
    b <- gl(2, 100)
    expect_equal(simpleSumFactors(y, block=b, min.mean=0), simpleSumFactors(as.matrix(y), block=b, min.mean=0))

    library(HDF5Array)
    x <- matrix(rgamma(ngenes*ncells, 2, 2), ngenes, ncells)
    z <- as(x, "HDF5Array")
    expect_equal(simpleSumFactors(z, min.mean=0), simpleSumFactors(x, min.mean=0))
    expect_equal(simpleSumFactors(z, block=b, min.mean=0), simpleSumFactors(x, block=b, min.mean=0))
})

set.seed(120005)
test_that("simpleSumFactors() behaves gracefully when encountering nonsensical size factors", {
    ncells <- 200
    ngenes <- 1000
    x <- matrix(rpois(ngenes*ncells, lambda=runif(100)), nrow=ngenes, ncol=ncells, byrow=TRUE)

    # Cleaning of size factors is properly triggered.
    x2 <- x
    x2[,1:20] <- 0 # large distance to 20th neighbor from any of these cells, so this will not be the reference.
    x2[1,1:20] <- 100
    expect_warning(out <- simpleSumFactors(x2, min.mean=0), "zero")
    expect_true(all(out > 0))

    x2 <- x
    x2[,1:21] <- 0 # one of these cells gets chosen as the reference.
    x2[1,1:21] <- 100
    expect_error(suppressWarnings(simpleSumFactors(x2, min.mean=0)), "infinite size factors")

    # Zero rescaling across-block handling.
    x2 <- x
    x2[,1:30] <- 0
    x2[1:100,1:30] <- rpois(3000, lambda=20)
    expect_error(out <- simpleSumFactors(x2, min.mean=0, block=rep(1:2, c(30, ncells-30))), "not strictly positive")

    # Filter threshold too strong.
    expect_error(simpleSumFactors(x, min.mean=10), "no genes")

    # No cells.
    expect_error(simpleSumFactors(x[,0]), "no neighbors")
})
