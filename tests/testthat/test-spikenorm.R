# Checking out the behaviour of the computeSpikeFactors function.
# library(testthat); library(scran); source("test-spikenorm.R")

set.seed(20003)
ncells <- 200
ngenes <- 1000

dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
is.spike <- rbinom(ngenes, 1, 0.7)==0L
dummy[is.spike,] <- matrix(rnbinom(sum(is.spike)*ncells, mu=20, size=5), ncol=ncells, nrow=sum(is.spike), byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
X <- splitSCEByAlt(X, ifelse(is.spike, "ERCC", "gene"))

test_that("computeSpikeFactors calculates spike-based size factors correctly", {
    out <- computeSpikeFactors(X, "ERCC")
    ref <- colSums(dummy[is.spike,])
    expect_equal(unname(sizeFactors(out)), ref/mean(ref))
})

test_that("computeSpikeFactors fails correctly on silly inputs", {
    expect_error(computeSpikeFactors(X, "ERR"))

    # Checking that it correctly returns nothing.
    out <- computeSpikeFactors(X[,0], "ERCC")
    expect_identical(unname(sizeFactors(out)), numeric(0))
})
