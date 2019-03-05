# This tests the multi-block normalization and variance calculations.
# require(scran); require(testthat); source("test-multiblock.R")

set.seed(20009)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
isSpike(X, "MySpikes") <- rbinom(ngenes, 1, 0.7)==0
sizeFactors(X, "MySpikes") <- runif(ncells)

block <- sample(5, ncells, replace=TRUE)

test_that("multiBlockNorm works properly", {
    X0 <- normalize(X)
    X2 <- multiBlockNorm(X, block)

    # Checking size factors.
    expect_identical(sizeFactors(X2), sizeFactors(X0))
    expect_false(!isTRUE(all.equal(sizeFactors(X2, "MySpike"), sizeFactors(X0, "MySpike"))))

    ref.en <- sizeFactors(X2)
    sp <- sizeFactors(X2, "MySpikes")
    for (i in unique(block)) {
        chosen <- block==i
        target <- mean(ref.en[chosen])
        expect_equal(target, mean(sp[chosen]))

        original <- sizeFactors(X, "MySpikes")[chosen]
        expect_equal(sp[chosen], original/mean(original)*target)
    }

    # Checking that log-counts are correctly computed.
    expect_equal(logcounts(X2)[!isSpike(X),], log2(t(t(counts(X)[!isSpike(X),])/sizeFactors(X2))+1))
    expect_equal(logcounts(X2)[isSpike(X),], log2(t(t(counts(X)[isSpike(X),])/sizeFactors(X2, "MySpikes"))+1))

    # Checking that library size substitution works.
    Y <- X
    sizeFactors(Y) <- NULL
    expect_warning(Y <- multiBlockNorm(Y, block), "using library sizes")

    Y2 <- X
    sizeFactors(Y2) <- scater::librarySizeFactors(Y2)
    Y2 <- multiBlockNorm(Y2, block)
    expect_equal(logcounts(Y), logcounts(Y2))
    expect_equal(sizeFactors(Y, "MySpikes"), sizeFactors(Y2, "MySpikes"))

    # Checking that it does nothing when there are no spike-ins or their size factors.
    Y3 <- clearSpikes(X)
    Y3 <- multiBlockNorm(Y3, block=block)
    X3 <- normalize(Y3)
    expect_equal(X3, Y3)
})

test_that("multiBlockVar works properly", {
    X2 <- multiBlockNorm(X, block)
    multi <- multiBlockVar(X2, block=block)

    expect_identical(colnames(multi$per.block), as.character(sort(unique(block))))
    expect_identical(rownames(multi), rownames(X))

    # Checking that the subsetted values were correct.
    for (b in unique(block)) {
        X3 <- X2[,block==b]
        fit <- trendVar(X3)
        dec <- decomposeVar(X3, fit)

        ref <- multi$per.block[[as.character(b)]]
        expect_equal(metadata(ref)$trend(1:5), fit$trend(1:5))
        metadata(ref)$trend <- NULL
        expect_equal(as.data.frame(dec), as.data.frame(ref))
    }

    # Checking that the combined variances are correct.
    per.block <- multi$per.block
    multi2 <- multi
    multi2$per.block <- NULL
    expect_identical(multi2, do.call(combineVar, as.list(per.block)))

    # Checking that it works happily with different settings.
    multi2 <- multiBlockVar(X2, block=block, trend.args=list(parametric=TRUE))
    expect_identical(nrow(multi), nrow(multi2))
    expect_identical(colnames(multi), colnames(multi2))

    multi2 <- multiBlockVar(X2, block=block, method="z")
    expect_identical(nrow(multi), nrow(multi2))
    expect_identical(colnames(multi), colnames(multi2))

    multi3 <- multiBlockVar(X2, block=block, make.tech.trend=TRUE)
    expect_identical(nrow(multi), nrow(multi3))
    expect_identical(colnames(multi), colnames(multi3))

    # Expecting a warning.
    Y <- normalize(X)
    expect_warning(multi <- multiBlockVar(Y, block=block), "centred")

    silly.block <- c(1,2,rep(3, ncol(X2)-2))
    expect_warning(multi <- multiBlockVar(X2, block=silly.block), "fewer than two cells")
    expect_error(multi <- multiBlockVar(X2, block=seq_len(ncol(X2))), "no block")
})
