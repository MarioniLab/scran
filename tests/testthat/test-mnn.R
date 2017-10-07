# Checks the application of the mnnCorrect function.
# require(scran); require(testthat); source("test-mnn.R")

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
