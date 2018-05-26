# This tests the functions related to mnnCorrect2.
# library(scran); library(testthat); source("test-mnn2.R")

expect_equal_besides_sign <- function(left, right, ...) {
    ratio <- left/right
    right2 <- sweep(right, 2, sign(ratio[1,]), FUN="*")
    expect_equal(left, right2, ...)
}

set.seed(1200001)
test_that("multi-sample PCA works as expected", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)

    # Checking that we get the same result independent of the multiples of the data. 
    ref <- scran:::.multi_pca(list(test1, test2), d=5)
    out <- scran:::.multi_pca(list(test1, cbind(test2, test2)), d=5)
    expect_equal(ref[[1]], out[[1]])
    expect_equal(rbind(ref[[2]], ref[[2]]), out[[2]])
    expect_identical(ncol(ref[[1]]), 5L)
    expect_identical(ncol(ref[[2]]), 5L)

    # Another check. 
    ref <- scran:::.multi_pca(list(test1, test2), d=3)
    out <- scran:::.multi_pca(list(cbind(test1, test1, test1), test2), d=3)
    expect_equal(ref[[2]], out[[2]])
    expect_equal(rbind(ref[[1]], ref[[1]], ref[[1]]), out[[1]])
    expect_identical(ncol(ref[[1]]), 3L)
    expect_identical(ncol(ref[[2]]), 3L)

    # Checking that the distances match up to the original values when we use full rank.
    out <- scran:::.multi_pca(list(test1, test2), d=10)
    
    rx <- seq_len(ncol(test1))
    cx <- ncol(test1) + seq_len(ncol(test2))
    ref.dist <- as.matrix(dist(rbind(t(test1), t(test2))))[rx, cx]
    out.dist <- as.matrix(dist(rbind(out[[1]], out[[2]])))[rx, cx]
    expect_equal(ref.dist, out.dist)
    
    out2 <- scran:::.multi_pca(list(test1, cbind(test2, test2)), d=10)
    out.dist2 <- as.matrix(dist(rbind(out2[[1]], out2[[2]])))[rx, cx]
    expect_equal(out.dist2, out.dist)

    # Checking that we get similar results with irlba.
    ref <- scran:::.multi_pca(list(test1, test2), d=2)
    out <- scran:::.multi_pca(list(test1, test2), d=2, approximate=TRUE)
    expect_equal_besides_sign(out[[1]], ref[[1]], tol=1e-7)
    expect_equal_besides_sign(out[[2]], ref[[2]], tol=1e-7)
})

set.seed(1200002)
test_that("averaging correction vectors works as expected", {
    test1 <- matrix(rnorm(1000), ncol=10)
    test2 <- matrix(rnorm(2000), ncol=10)
    mnn1 <- sample(nrow(test1), 250, replace=TRUE)
    mnn2 <- sample(nrow(test1), 250, replace=TRUE)

    # Slow reference calculation.
    correct <- test1[mnn1,] - test2[mnn2,]
    by.mnn <- split(seq_along(mnn2), mnn2)
    collected <- vector("list", length(by.mnn))
    for (idx in seq_along(by.mnn)) {
        collected[[idx]] <- colMeans(correct[by.mnn[[idx]],,drop=FALSE])
    }
    ref <- do.call(rbind, collected)
    rownames(ref) <- names(by.mnn)

    # Comparing the implementation in scran.
    out <- scran:::.average_correction(test1, mnn1, test2, mnn2)  
    expect_equal(out$averaged, ref)
    expect_identical(out$second, sort(unique(mnn2)))
    expect_identical(as.integer(rownames(out$averaged)), out$second)
})

set.seed(1200003)
test_that("centering along a batch vector works correctly", {
    test <- matrix(rnorm(1000), ncol=10)
    batch <- rnorm(10) 
    centered <- scran:::.center_along_batch_vector(test, batch)
    new.locations <- centered %*% batch
    expect_true(mad(new.locations) < 1e-8)
})

