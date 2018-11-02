# Checks the application of the mnnCorrect function.
# require(scran); require(testthat); source("test-mnn.R")

set.seed(10001)
test_that("Cosine normalization is correct", {
    X <- matrix(rnorm(10000), ncol=100)
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, scran:::cosineNorm(X))

    X <- matrix(rpois(20000, lambda=5), ncol=1000)
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, scran:::cosineNorm(X))
})

set.seed(10000)
test_that("Mutual NN detection is correct", {
    # Reference NNs.
    library(Matrix)
    REF <- function(d1, d2, k1, k2) {
        n1 <- nrow(d1)
        n2 <- nrow(d2)
        n.total <- n1 + n2
   
        W21 <- BiocNeighbors::queryKNN(d2, query=d1, k=k1)
        W12 <- BiocNeighbors::queryKNN(d1, query=d2, k=k2)
        W <- sparseMatrix(i=c(rep(seq_len(n1), k1), rep(n1 + seq_len(n2), k2)),
                          j=c(n1 + W21$index, W12$index),
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
    comparator(REF(A, B, 10, 10), scran:::find.mutual.nn(A, B, 10, 10))
    comparator(REF(A, B, 5, 20), scran:::find.mutual.nn(A, B, 5, 20))
    comparator(REF(A, B, 20, 5), scran:::find.mutual.nn(A, B, 20, 5))

    A <- matrix(rpois(25000, lambda=20), ncol=100)
    B <- matrix(rpois(15000, lambda=50), ncol=100)
    comparator(REF(A, B, 10, 10), scran:::find.mutual.nn(A, B, 10, 10))
    comparator(REF(A, B, 5, 20), scran:::find.mutual.nn(A, B, 5, 20))
    comparator(REF(A, B, 20, 5), scran:::find.mutual.nn(A, B, 20, 5))
})

set.seed(10002)
test_that("Biological subspace is correctly re-projected", {
    A <- matrix(rnorm(10000), ncol=50)
    subset <- sample(nrow(A), nrow(A)/2)   
    ref <- scran:::get.bio.span(A[subset,], 3) 
    proj <- scran:::get.bio.span(A, 3, subset.row=subset) 
    expect_equal(ref, proj[subset,])

    B <- rbind(A, A)
    first.half <- seq_len(nrow(A))
    ref <- scran:::get.bio.span(A, 3)
    proj <- scran:::get.bio.span(B, 3, subset.row=first.half) 
    expect_equal(ref, proj[first.half,])
    expect_equal(proj[first.half,], proj[-first.half,])
})

set.seed(10003)
test_that("Batch vectors are correctly calculated", {
    data1 <- matrix(rnorm(10000, sd=0.1), ncol=25)
    data2 <- matrix(rnorm(25000, sd=0.1), ncol=25)
    mnn1 <- 1:10
    mnn2 <- 30:21

    # Constructing a reference function.
    REF <- function(data1, data2, mnn1, mnn2, s2) {
        d <- as.matrix(dist(data2))
        w <- exp(-d^2/s2)

        mnn.dens <- rowSums(w[,unique(mnn2)])
        N <- tabulate(mnn2, nbins=nrow(data2))
        kernel <- t(w/(N*mnn.dens))[,mnn2]

        kernel <- kernel/rowSums(kernel)
        vect <- data1[mnn1,] - data2[mnn2,]
        out <- kernel %*% vect
        dimnames(out) <- NULL
        return(out)
    }

    # Vanilla check
    s2 <- 0.1
    xx <- scran:::compute.correction.vectors(data1, data2, mnn1, mnn2, t(data2), s2)
    ref <- REF(data1, data2, mnn1, mnn2, s2)
    expect_equal(xx, ref)

    # Check with cells in multiple MNN pairs.
    alt.mnn1 <- c(11, 12, 13, mnn1)
    alt.mnn2 <- c(30, 30, 30, mnn2)
    xx <- scran:::compute.correction.vectors(data1, data2, alt.mnn1, alt.mnn2, t(data2), s2)
    ref <- REF(data1, data2, alt.mnn1, alt.mnn2, s2)
    expect_equal(xx, ref)

    # Check with more MNN pairs involved.
    alt.mnn1 <- 1:200
    alt.mnn2 <- 500:301
    xx <- scran:::compute.correction.vectors(data1, data2, alt.mnn1, alt.mnn2, t(data2), s2)
    ref <- REF(data1, data2, alt.mnn1, alt.mnn2, s2)
    expect_equal(xx, ref)

    # Check with a different bandwidth. 
    s2 <- 0.5
    xx <- scran:::compute.correction.vectors(data1, data2, mnn1, mnn2, t(data2), s2)
    ref <- REF(data1, data2, mnn1, mnn2, s2)
    expect_equal(xx, ref)
})

set.seed(100032)
test_that("Variance shift adjustment is correctly performed", {
    data1 <- matrix(rnorm(10000, sd=0.1), ncol=25)
    data2 <- matrix(rnorm(25000, sd=0.1), ncol=25)
    corvect <- matrix(runif(length(data2)), nrow=nrow(data2))

    # Constructing a reference function.
    REF <- function(data1, data2, cell.vect, sigma) {
        scaling <- numeric(nrow(cell.vect))
        for (cell in seq_along(scaling)) {
            # For each cell, projecting both data sets onto the normalized correction vector for that cell.
            cur.cor.vect <- cell.vect[cell,]
            l2norm <- sqrt(sum(cur.cor.vect^2))
            cur.cor.vect <- cur.cor.vect/l2norm
            coords2 <- data2 %*% cur.cor.vect
            coords1 <- data1 %*% cur.cor.vect
    
            # Also getting the distance from the correction vector. 
            dist2 <- data2[cell,] - t(data2)
            dist2 <- dist2 - outer(cur.cor.vect, as.numeric(crossprod(dist2, cur.cor.vect)))
            dist2 <- colSums(dist2^2)
            weight2 <- exp(-dist2/sigma)
    
            dist1 <- data2[cell,] - t(data1)
            dist1 <- dist1 - outer(cur.cor.vect, as.numeric(crossprod(dist1, cur.cor.vect)))
            dist1 <- colSums(dist1^2)
            weight1 <- exp(-dist1/sigma)
    
            # Computing the weighted cumulative probability, for quantile-quantile mapping.
            rank2 <- rank(coords2, ties.method="first")
            prob2 <- sum(weight2[rank2 <= rank2[cell]])/sum(weight2)
            ord1 <- order(coords1)
            ecdf1 <- cumsum(weight1[ord1])/sum(weight1)
            
            # Adjusting the length of the correction vector so that the correction will match the quantiles.
            quan1 <- coords1[ord1[min(which(ecdf1 >= prob2))]]
            quan2 <- coords2[cell]
            scaling[cell] <- (quan1 - quan2)/l2norm
        }
        return(scaling)
    }

    ref <- REF(data1, data2, corvect, 1)
    test <- .Call(scran:::cxx_adjust_shift_variance, t(data1), t(data2), corvect, 1)
    expect_equal(ref, test)

    ref <- REF(data1, data2, corvect, 0.1)
    test <- .Call(scran:::cxx_adjust_shift_variance, t(data1), t(data2), corvect, 0.1)
    expect_equal(ref, test)
})

set.seed(10004)
test_that("mnnCorrect behaves consistently with subsetting", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)

    keep <- 1:5 
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)    

    # Without cosine normalization of the output.
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], cos.norm.out=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, cos.norm.out=FALSE)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)    

    # Without cosine normalization of the input.
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], cos.norm.in=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, cos.norm.in=FALSE)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)    

    # Without any cosine normalization at all.
    keep <- 6:10
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], cos.norm.in=FALSE, cos.norm.out=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, cos.norm.in=FALSE, cos.norm.out=FALSE)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)    

    # With SVDs.
    keep <- 2:7
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], svd.dim=2)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, svd.dim=2)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)   

    # Without the variance adjustment.
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], var.adj=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, var.adj=FALSE)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)   

    # Asking for angles.
    keep <- 8:3
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], compute.angle=TRUE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, compute.angle=TRUE)
    out$corrected <- lapply(out$corrected, "[", i=keep,)
    expect_equal(ref, out)   

    # Duplicated genes should have no effect.
    out <- mnnCorrect(rbind(alpha, alpha), rbind(bravo, bravo), rbind(charlie, charlie), 
                      subset.row=1:nrow(alpha), svd.dim=2)
    ref1 <- lapply(out$corrected, "[", i=1:nrow(alpha),)
    ref2 <- lapply(out$corrected, "[", i=nrow(alpha)+1:nrow(alpha),)
    expect_equal(ref1, ref2) 
    
    # Checking that the order has no effect.
    new.order <- c(2, 3, 1)
    out <- mnnCorrect(alpha, bravo, charlie, order=new.order)
    ref <- mnnCorrect(bravo, charlie, alpha)
    expect_equal(out$corrected[new.order], ref$corrected)
})

set.seed(10004)
library(Matrix)
test_that("mnnCorrect behaves properly with sparse matrices", {
    alpha <- rsparsematrix(20, 100, density=0.5)
    bravo <- rsparsematrix(20, 100, density=0.5)
    charlie <- rsparsematrix(20, 100, density=0.5)

    out <- mnnCorrect(alpha, bravo, charlie)
    ref <- mnnCorrect(as.matrix(alpha), as.matrix(bravo), as.matrix(charlie))
    expect_equivalent(ref$corrected, lapply(out$corrected, as.matrix))
    expect_identical(out$pairs, ref$pairs)
})

