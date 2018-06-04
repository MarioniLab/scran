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
isSpike(X, "MySpike") <- is.spike

test_that("computeSpikeFactors calculates spike-based size factors correctly", {
    out <- computeSpikeFactors(X)
    ref <- colSums(dummy[is.spike,])
    expect_equal(unname(sizeFactors(out)), ref/mean(ref))
    expect_equal(sizeFactors(out), sizeFactors(out, type="MySpike"))
})

test_that("computeSpikeFactors works with multiple spike-in sets", {
    subset <- split(which(is.spike), rep(1:2, length.out=sum(is.spike)))
    isSpike(X, "MySpike") <- subset[[1]]
    isSpike(X, "SecondSpike") <- subset[[2]]
    out <- computeSpikeFactors(X)

    out.sub <- computeSpikeFactors(X, type="MySpike") # Sanity check, to make sure that it's calculating it differently for each spike-in.
    subref <- colSums(dummy[subset[[1]],])
    expect_equal(unname(sizeFactors(out.sub)), subref/mean(subref))
    expect_equal(sizeFactors(out.sub), sizeFactors(out.sub, type="MySpike"))
    expect_identical(sizeFactors(out.sub, type="SecondSpike"), NULL)
    
    out2 <- computeSpikeFactors(X)
    expect_equal(sizeFactors(out), sizeFactors(out2))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))
    
    out2 <- computeSpikeFactors(X, type=c("MySpike", "SecondSpike"))
    expect_equal(sizeFactors(out), sizeFactors(out2))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))
})

test_that("computeSpikeFactors responds correctly to general.use", {
    sizeFactors(X) <- 1
    out <- computeSpikeFactors(X, general.use=FALSE)
    expect_equal(unname(sizeFactors(out)), rep(1, ncells))
    ref <- colSums(dummy[is.spike,])
    expect_equal(unname(sizeFactors(out, type="MySpike")), ref/mean(ref))
})

test_that("computeSpikeFactors fails correctly on silly inputs", {
    expect_error(out <- computeSpikeFactors(X[0,]), "no spike-in transcripts present in 'x'")

    alt.X <- X
    counts(alt.X)[] <- 0L
    expect_warning(out <- computeSpikeFactors(alt.X), "zero spike-in counts during spike-in normalization")
    expect_identical(unname(sizeFactors(out)), rep(NaN, ncol(out)))

    # Checking that it correctly returns nothing.
    out <- computeSpikeFactors(X[,0])
    expect_identical(unname(sizeFactors(out)), numeric(0))
})
