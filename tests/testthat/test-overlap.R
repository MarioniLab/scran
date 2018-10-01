# This checks the overlapExprs function.
# require(scran); require(testthat); source("test-overlap.R")

set.seed(11000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
grouping <- as.character(sample(3, Ncells, replace=TRUE))

test_that("overlapExprs behaves consistently with SingleCellExperiment objects", {
    X2 <- SingleCellExperiment(list(counts=X))
    sizeFactors(X2) <- colSums(X)
    X2 <- normalize(X2)

    grouping <- rep(1:4, each=25)
    expect_equal(overlapExprs(exprs(X2), grouping), overlapExprs(X2, grouping)) 
    expect_equal(overlapExprs(counts(X2), grouping), overlapExprs(X2, grouping, assay.type="counts")) 

    isSpike(X2, "MySpike") <- rbinom(Ngenes, 1, 0.6)==0L
    expect_equal(overlapExprs(exprs(X2), grouping, subset.row=!isSpike(X2)), overlapExprs(X2, grouping)) 
    expect_equal(overlapExprs(exprs(X2), grouping), overlapExprs(X2, grouping, get.spikes=TRUE)) 
})

test_that("overlapExprs fails correctly on silly examples", {
    #  No genes.
    grouping <- rep(1:4, each=25)
    out <- overlapExprs(X, grouping, subset.row=integer(0))
    expect_identical(names(out), as.character(1:4))
    expect_identical(unname(sapply(out, nrow)), integer(length(out))) 
    out2 <- overlapExprs(X[0,], grouping)
    expect_identical(out, out2)
    
    # No cells, or mismatched numbers of cells.
    expect_error(overlapExprs(X[,0], grouping[0]), "need at least two")
    expect_error(overlapExprs(X, grouping[-1]), "length of 'clusters' does not equal 'ncol(x)'", fixed=TRUE)
})

