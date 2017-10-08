# Checks the application of the mnnCorrect function.
# require(scran); require(testthat); source("test-mnn.R")

set.seed(10001)
test_that("Cosine normalization is correct", {
    X <- matrix(rnorm(10000), ncol=100)
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, scran:::cosine.norm(X))

    X <- matrix(rpois(20000, lambda=5), ncol=1000)
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, scran:::cosine.norm(X))
})

set.seed(10000)
test_that("Mutual NN detection is correct", {
    # Reference NNs.
    library(FNN)
    library(Matrix)
    REF <- function(d1, d2, k1, k2) {
        n1 <- nrow(d1)
        n2 <- nrow(d2)
        n.total <- n1 + n2
   
        W21 <- get.knnx(d2, query=d1, k=k1)
        W12 <- get.knnx(d1, query=d2, k=k2)
        W <- sparseMatrix(i=c(rep(seq_len(n1), k1), rep(n1 + seq_len(n2), k2)),
                          j=c(n1 + W21$nn.index, W12$nn.index),
                          x=rep(1, n1*k1 + n2*k2), dims=c(n.total, n.total))

        W <- W * t(W) # elementwise multiplication to keep mutual nns only
        A <- which(W>0, arr.ind=TRUE) # row/col indices of mutual NNs

        A1 <- A[,1]
        A1 <- A1[A1 <= n1]
        A2 <- A[,2] - n1
        A2 <- A2[A2 > 0]
        return(list(first=A1, second=A2))
    }
    
    # Check that the values are identical.
    comparator <- function(x, y) {
        ox <- order(x$first, x$second)
        oy <- order(y$first, y$second)
        expect_identical(x$first[ox], y$first[oy])
        expect_identical(x$second[ox], y$second[oy]) 
    }
    
    # Compare to actual run.
    A <- matrix(rnorm(10000), ncol=50)
    B <- matrix(rnorm(20000), ncol=50)
    comparator(REF(A, B, 10, 10), scran:::find.mutual.nn(A, B, 10, 10, SerialParam()))
    comparator(REF(A, B, 10, 20), scran:::find.mutual.nn(A, B, 10, 20, SerialParam()))
    comparator(REF(A, B, 5, 20), scran:::find.mutual.nn(A, B, 5, 20, SerialParam()))
    comparator(REF(A, B, 20, 5), scran:::find.mutual.nn(A, B, 20, 5, SerialParam()))

    A <- matrix(rpois(25000, lambda=20), ncol=100)
    B <- matrix(rpois(15000, lambda=50), ncol=100)
    comparator(REF(A, B, 10, 10), scran:::find.mutual.nn(A, B, 10, 10, SerialParam()))
    comparator(REF(A, B, 10, 20), scran:::find.mutual.nn(A, B, 10, 20, SerialParam()))
    comparator(REF(A, B, 5, 20), scran:::find.mutual.nn(A, B, 5, 20, SerialParam()))
    comparator(REF(A, B, 20, 5), scran:::find.mutual.nn(A, B, 20, 5, SerialParam()))
})

set.seed(10002)
test_that("Smoothing kernel construction is correct", {
    data <- matrix(rnorm(10000, sd=0.1), ncol=25)
    
    # No kernel.
    expect_identical(matrix(1, nrow(data), nrow(data)), scran:::construct.smoothing.kernel(data, sigma=NA))

    # Full gaussian kernel.
    d <- as.matrix(dist(data))
    s <- 0.1
    out <- scran:::construct.smoothing.kernel(data, sigma=s)
    expect_equal(out, exp(-d^2/s))
    s <- 0.5
    out <- scran:::construct.smoothing.kernel(data, sigma=s)
    expect_equal(out, exp(-d^2/s))

    # Truncated kernel.
    REF <- function(data, sigma, kk, mnn.set) {
        mnn.set <- unique(mnn.set)
        nn.out <- get.knnx(data[mnn.set,,drop=FALSE], query=data, k=kk)
        result <- matrix(0, nrow(data), nrow(data))
        for (i in seq_len(nrow(data))) {
            result[i,mnn.set[nn.out$nn.index[i,]]] <- exp(-nn.out$nn.dist[i,]^2/sigma)
        }
        return(result)
    }

    mnn.set <- sample(nrow(data), 200)
    xx <- as.matrix(scran:::construct.smoothing.kernel(data, sigma=0.1, kk=100, mnn.set=mnn.set, exact=FALSE, BPPARAM=SerialParam()))
    dimnames(xx) <- NULL
    expect_equal(xx, REF(data, sigma=0.1, kk=100, mnn.set=mnn.set))
    
    mnn.set <- c(1:20, 10:50)
    xx <- as.matrix(scran:::construct.smoothing.kernel(data, sigma=0.1, kk=50, mnn.set=mnn.set, exact=FALSE, BPPARAM=SerialParam()))
    dimnames(xx) <- NULL
    expect_equal(xx, REF(data, sigma=0.1, kk=50, mnn.set=mnn.set))
  
    mnn.set <- c(100:20, 10:50)
    xx <- as.matrix(scran:::construct.smoothing.kernel(data, sigma=0.1, kk=25, mnn.set=mnn.set, exact=FALSE, BPPARAM=SerialParam()))
    dimnames(xx) <- NULL
    expect_equal(xx, REF(data, sigma=0.1, kk=25, mnn.set=mnn.set)) 
})

set.seed(10003)
test_that("Batch vectors are correctly calculated", {
    data1 <- matrix(rnorm(10000, sd=0.1), ncol=25)
    data2 <- matrix(rnorm(25000, sd=0.1), ncol=25)
    mnn1 <- 1:10
    mnn2 <- 30:21

    # Using a full kernel for an easier check.
    kernel <- matrix(1, nrow(data2), nrow(data2))
    xx <- scran:::compute.correction.vectors(data1, data2, mnn1, mnn2, kernel)
    ref <- data1[mnn1,] - data2[mnn2,]
    vec <- colMeans(ref)
    expect_equal(xx, matrix(vec, nrow(xx), ncol(xx), byrow=TRUE))

    # Adding redundancies shouldn't change the result.
    xx <- scran:::compute.correction.vectors(data1, data2, c(1,1,1,1,mnn1), c(30,30,30,30,mnn2), kernel)
    expect_equal(xx, matrix(vec, nrow(xx), ncol(xx), byrow=TRUE))

    # Scaling the kernel shouldn't change the result.
    xx <- scran:::compute.correction.vectors(data1, data2, mnn1, mnn2, kernel*2)
    expect_equal(xx, matrix(vec, nrow(xx), ncol(xx), byrow=TRUE))

    # Trying out a somewhat less trivial kernel.
    data1 <- matrix(rnorm(10000, sd=0.1), ncol=25)
    data2 <- matrix(rnorm(25000, sd=0.1), ncol=25)
    mnn1 <- 10:40
    mnn2 <- 30:60
    ref <- data1[mnn1,] - data2[mnn2,]

    kernel <- matrix(runif(nrow(data2)^2), nrow(data2), nrow(data2))
    norm.kernel <- kernel[,mnn2]
    norm.kernel <- norm.kernel/rowSums(norm.kernel)
    vec2 <- norm.kernel %*% ref
    expect_equal(vec2, scran:::compute.correction.vectors(data1, data2, mnn1, mnn2, kernel))
    expect_equal(vec2, scran:::compute.correction.vectors(data1, data2, c(10, 10, 10, mnn1), c(30, 30, 30, mnn2), kernel))
})
