# This tests the sandbag function, by running it and checking that the selecter markers make sense.
# require(scran); require(testthat); source("test-sandbag.R")

happycheck <- function(X1, X2, X3, pairings, frac) {
    is.okay <- logical(nrow(pairings))
    thresh1 <- ceiling(ncol(X1) * frac)
    thresh2 <- ceiling(ncol(X2) * frac)
    thresh3 <- ceiling(ncol(X3) * frac)
    for (p in seq_along(is.okay)) {
         diff1 <- X1[pairings$first[p],] - X1[pairings$second[p],]
         diff2 <- X2[pairings$first[p],] - X2[pairings$second[p],]
         diff3 <- X3[pairings$first[p],] - X3[pairings$second[p],]
         
         u1 <- sum(diff1 > 0)
         u2 <- sum(diff2 > 0)
         u3 <- sum(diff3 > 0)
         d1 <- sum(diff1 < 0)
         d2 <- sum(diff2 < 0)
         d3 <- sum(diff3 < 0)

         is.okay[p] <- (u1 >= thresh1 && d2 >= thresh2 && d3 >= thresh3) || 
                       (d1 >= thresh1 && u2 >= thresh2 && u3 >= thresh3)
    }
    return(is.okay)
}

####################################################################################################

set.seed(100)

Ngenes <- 100
phases <- sample(3, 100, replace=TRUE)
is.G1 <- phases==1L
is.G2M <- phases==2L
is.S <- phases==3L

frac <- 0.5
X <- matrix(rpois(Ngenes*length(phases), lambda=10), nrow=Ngenes)
rownames(X) <- paste0("X", seq_len(Ngenes))
cur.classes <- list(G1=is.G1, S=is.S, G2M=is.G2M)
out <- sandbag(X, cur.classes, fraction=frac)

XG1 <- X[,is.G1,drop=FALSE]
XS <- X[,is.S,drop=FALSE]
XG2M <- X[,is.G2M,drop=FALSE]

expect_true(all(happycheck(XG1, XS, XG2M, out$G1, frac)))
expect_true(all(happycheck(XG2M, XS, XG1, out$G2M, frac)))
expect_true(all(happycheck(XS, XG1, XG2M, out$S, frac)))

# Checking silly inputs.

out <- sandbag(X[0,], cur.classes, fraction=frac)
expect_identical(out$G1, data.frame(first=character(0), second=character(0), stringsAsFactors=FALSE))
expect_identical(out$G2M, data.frame(first=character(0), second=character(0), stringsAsFactors=FALSE))
expect_identical(out$S, data.frame(first=character(0), second=character(0), stringsAsFactors=FALSE))
expect_error(sandbag(X, list(G1=integer(0), S=is.S, G2M=is.G2M), fraction=frac), "each class must have at least one cell")
expect_error(sandbag(X, unname(cur.classes), fraction=frac), "names")

is.G1 <- 1
is.G2M <- 2
is.S <- 3
out <- sandbag(X, list(G1=is.G1, S=is.S, G2M=is.G2M), fraction=frac)
XG1 <- X[,is.G1,drop=FALSE]
XS <- X[,is.S,drop=FALSE]
XG2M <- X[,is.G2M,drop=FALSE]
expect_true(all(happycheck(XG1, XS, XG2M, out$G1, frac)))
expect_true(all(happycheck(XG2M, XS, XG1, out$G2M, frac)))
expect_true(all(happycheck(XS, XG1, XG2M, out$S, frac)))

# Testing for a SCESet, without spike-ins.

set.seed(200)
test_that("sandbag works correctly with SingleCellExperiment objects", {
    Ngenes <- 100
    phases <- sample(3, 100, replace=TRUE)
    is.G1 <- phases==1L
    is.G2M <- phases==2L
    is.S <- phases==3L
    
    frac <- 0.5
    X <- matrix(rpois(Ngenes*length(phases), lambda=10), nrow=Ngenes)
    rownames(X) <- paste0("X", seq_len(Ngenes))
    X <- SingleCellExperiment(list(counts=X))
    
    cur.classes <- list(G1=is.G1, S=is.S, G2M=is.G2M)
    out <- sandbag(X, cur.classes, fraction=frac)
    expect_identical(sandbag(counts(X), cur.classes, fraction=frac), out)
    
    XG1 <- counts(X[,is.G1,drop=FALSE])
    XS <- counts(X[,is.S,drop=FALSE])
    XG2M <- counts(X[,is.G2M,drop=FALSE])
    
    expect_true(all(happycheck(XG1, XS, XG2M, out$G1, frac)))
    expect_true(all(happycheck(XG2M, XS, XG1, out$G2M, frac)))
    expect_true(all(happycheck(XS, XG1, XG2M, out$S, frac)))
})
