# This tests the functions related to multiBatchPCA.
# library(scran); library(testthat); source("test-multipca.R")

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
    expect_equal_besides_sign(ref[[1]], out[[1]])
    expect_equal_besides_sign(rbind(ref[[2]], ref[[2]]), out[[2]])
    expect_identical(ncol(ref[[1]]), 5L)
    expect_identical(ncol(ref[[2]]), 5L)

    # Another check. 
    ref <- scran:::.multi_pca(list(test1, test2), d=3)
    out <- scran:::.multi_pca(list(cbind(test1, test1, test1), test2), d=3)
    expect_equal_besides_sign(ref[[2]], out[[2]])
    expect_equal_besides_sign(rbind(ref[[1]], ref[[1]], ref[[1]]), out[[1]])
    expect_identical(ncol(ref[[1]]), 3L)
    expect_identical(ncol(ref[[2]]), 3L)

    # Checking with equal numbers of cells - should be equivalent to cbind'd PCA. 
    test3 <- matrix(rnorm(1000), nrow=10)
    out <- scran:::.multi_pca(list(test1, test3), d=4)
    ref <- prcomp(t(cbind(test1, test3)), rank.=4)
    expect_equal_besides_sign(rbind(out[[1]], out[[2]]), unname(ref$x))
    
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
})

set.seed(12000011)
test_that("fast SVD works as expected", {
    ngenes <- 50
    ncells1 <- 200
    test1 <- matrix(rnorm(ncells1 * ngenes), nrow=ngenes)
    ncells2 <- 400
    test2 <- matrix(rnorm(ncells2 * ngenes), nrow=ngenes)

    # Matches up well with reference irlba.
    t1 <- t(test1)
    t2 <- t(test2)
    com <- rbind(t1, t2)
    expect_true(nrow(com) > ncol(com))

    ref <- svd(com, nu=0, nv=5) 
    out <- scran:::.fast_svd(list(t1, t2), nv=5) 
    expect_equal(out$d, ref$d[1:5])
    expect_equal_besides_sign(out$v, ref$v)

    # Checking that we get similar results when nrows < ncols. 
    com2 <- rbind(t1[1:10,], t2[1:5,])
    expect_true(nrow(com2) < ncol(com2))
    ref <- svd(com2, nu=0, nv=5) 
    out <- scran:::.fast_svd(list(t1[1:10,], t2[1:5,]), nv=5) 
    expect_equal(out$d, ref$d[1:5])
    expect_equal_besides_sign(out$v, ref$v)

    # Checking that we get similar results from the whole suite.
    ref <- scran:::.multi_pca(list(test1, test2), d=2, use.crossprod=TRUE)
    out <- scran:::.multi_pca(list(test1, test2), d=2, use.crossprod=FALSE)
    expect_equal_besides_sign(out[[1]], ref[[1]])
    expect_equal_besides_sign(out[[2]], ref[[2]])
})

set.seed(12000012)
test_that("irlba works as expected", {
    ngenes <- 50
    ncells1 <- 250
    test1 <- matrix(rnorm(ncells1 * ngenes), nrow=ngenes)
    ncells2 <- 400
    test2 <- matrix(rnorm(ncells2 * ngenes), nrow=ngenes)

    # Matches up well with reference irlba.
    t1 <- t(test1)
    t2 <- t(test2)
    com <- rbind(t1, t2)

    set.seed(10)
    ref <- irlba::irlba(com, nu=0, nv=5) 
    set.seed(10)
    out <- scran:::.fast_svd(list(t1, t2), approximate=TRUE, nv=5) 
    expect_equal(out$d, ref$d)
    expect_equal_besides_sign(out$v, ref$v, tol=1e-4)

    set.seed(100)
    ref <- irlba::irlba(com, nu=0, nv=10) 
    set.seed(100)
    out <- scran:::.fast_svd(list(t1, t2), nv=10, approximate=TRUE) 
    expect_equal(out$d, ref$d)
    expect_equal_besides_sign(out$v, ref$v, tol=1e-4)

    # Forcing accuracy by increasing tolerance (also testing irlba.args!)
    expect_error(expect_equal_besides_sign(out$v, ref$v, tol=1e-8), "not equal")

    set.seed(1000)
    ref <- irlba::irlba(com, nu=0, nv=10, tol=1e-16) 
    set.seed(1000)
    out <- scran:::.fast_svd(list(t1, t2), nv=10, irlba.args=list(tol=1e-16), approximate=TRUE)
    expect_equal(out$d, ref$d)
    expect_equal_besides_sign(out$v, ref$v, tol=1e-8)

    # Checking that we get similar results when nrows < ncols. 
    set.seed(20)
    ref <- irlba::irlba(rbind(t1[1:10,], t2[1:5,]), nu=0, nv=5) 
    set.seed(20)
    out <- scran:::.fast_svd(list(t1[1:10,], t2[1:5,]), nv=5, approximate=TRUE)
    expect_equal(out$d, ref$d)
    expect_equal_besides_sign(out$v, ref$v, tol=1e-4)

    # Checking that we get similar results from the whole suite.
    set.seed(999)
    ref <- scran:::.multi_pca(list(test1, test2), d=2, approximate=TRUE, use.crossprod=TRUE)
    out <- scran:::.multi_pca(list(test1, test2), d=2, approximate=TRUE, use.crossprod=FALSE)
    expect_equal_besides_sign(out[[1]], ref[[1]], tol=1e-4)
    expect_equal_besides_sign(out[[2]], ref[[2]], tol=1e-4)

    ref <- scran:::.multi_pca(list(test1, test2), d=3, approximate=FALSE, use.crossprod=TRUE)
    out <- scran:::.multi_pca(list(test1, test2), d=3, approximate=TRUE, use.crossprod=TRUE)
    expect_equal_besides_sign(out[[1]], ref[[1]], tol=1e-4)
    expect_equal_besides_sign(out[[2]], ref[[2]], tol=1e-4)

    ref <- scran:::.multi_pca(list(test1, test2), d=4, approximate=FALSE, use.crossprod=FALSE)
    out <- scran:::.multi_pca(list(test1, test2), d=4, approximate=TRUE, use.crossprod=FALSE)
    expect_equal_besides_sign(out[[1]], ref[[1]], tol=1e-4)
    expect_equal_besides_sign(out[[2]], ref[[2]], tol=1e-4)
})
