# This tests the functionality of modelGeneVar.
# library(testthat); library(scran); source("test-model-var.R")

set.seed(20001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

library(scater)
dummy <- normalizeCounts(dummy)

test_that("modelGeneVar works correctly without blocking", {
    out <- modelGeneVar(dummy)
    expect_equal(out$mean, rowMeans(dummy))
    expect_equal(unname(out$total), DelayedMatrixStats::rowVars(dummy))
    expect_equal(out$tech, metadata(out)$trend(out$mean))
    expect_equal(out$bio, out$total-out$tech)
    expect_equal(order(out$p.value), order(out$bio/out$tech, decreasing=TRUE))
})

test_that("modelGeneVar works correctly with blocking, no weighting", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneVar(dummy, block=block)

    accumulated.mean <- accumulated.total <- accumulated.tech <- 0
    for (i in unique(block)) {
        ref <- modelGeneVar(dummy[,i==block])
        subout <- out$per.block[[i]]

        expect_identical(ref$mean, subout$mean)
        expect_identical(ref$total, subout$total)
        expect_identical(ref$tech, subout$tech)
        expect_identical(ref$bio, subout$bio)
        expect_identical(ref$p.value, subout$p.value)

        accumulated.mean <- accumulated.mean + ref$mean
        accumulated.total <- accumulated.total + ref$total
        accumulated.tech <- accumulated.tech + ref$tech
    }

    # Check combining statistics works correctly.
    n <- length(unique(block))
    expect_equal(out$mean, accumulated.mean/n)
    expect_equal(out$total, accumulated.total/n)
    expect_equal(out$tech, accumulated.tech/n)
    expect_equal(out$bio, out$total - out$tech)

    all.p <- lapply(out$per.block, "[[", i="p.value")
    expect_equal(out$p.value, do.call(combinePValues, all.p))

    # Responds to choice of method. 
    out2 <- modelGeneVar(dummy, block=block, method="z")
    all.p <- lapply(out2$per.block, "[[", i="p.value")
    expect_equal(out2$p.value, do.call(combinePValues, c(all.p, list(method='z'))))
})

test_that("modelGeneVar works correctly with blocking and weighting", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneVar(dummy, block=block, equiweight=FALSE)

    accumulated.mean <- accumulated.total <- accumulated.tech <- 0
    for (i in unique(block)) {
        ref <- modelGeneVar(dummy[,i==block])
        subout <- out$per.block[[i]]

        expect_identical(ref$mean, subout$mean)
        expect_identical(ref$total, subout$total)
        expect_identical(ref$tech, subout$tech)
        expect_identical(ref$bio, subout$bio)
        expect_identical(ref$p.value, subout$p.value)

        n <- sum(i==block)
        accumulated.mean <- accumulated.mean + ref$mean * n
        accumulated.total <- accumulated.total + ref$total * n
        accumulated.tech <- accumulated.tech + ref$tech * n
    }

    # Check combining statistics works correctly.
    n <- length(block)
    expect_equal(out$mean, accumulated.mean/n)
    expect_equal(out$total, accumulated.total/n)
    expect_equal(out$tech, accumulated.tech/n)
    expect_equal(out$bio, out$total - out$tech)

    all.p <- lapply(out$per.block, "[[", i="p.value")
    expect_equal(out$p.value, do.call(combinePValues, all.p))

    # Responds to choice of method with weighting.
    out2 <- modelGeneVar(dummy, block=block, method="z", equiweight=FALSE)
    all.p <- lapply(out2$per.block, "[[", i="p.value")
    w <- countMatches(names(all.p), block)
    expect_equal(out2$p.value, do.call(combinePValues, c(all.p, list(method='z', weights=w))))
})

test_that("modelGeneVar works with subsetting options", {
    chosen <- sample(ngenes, ngenes/2)
    out <- modelGeneVar(dummy, subset.row=chosen)
    ref <- modelGeneVar(dummy[chosen,])
    expect_equal(out, ref)

    # Subsetting by fit works.
    out2 <- modelGeneVar(dummy, subset.fit=chosen)
    expect_identical(rownames(out2), rownames(dummy))
    expect_equal(out2[chosen,1:5], ref[,1:5])

    # Zero-length subsetting.
    empty <- modelGeneVar(dummy, subset.row=integer(0), subset.fit=chosen)
    expect_identical(nrow(empty), 0L)

    expect_error(modelGeneVar(dummy, subset.fit=integer(0)), "need at least 2 points")
})

test_that("modelGeneVar works with design matrices", {
    # Testing with a trivial design matrix.
    design0 <- matrix(1, ncol(dummy), 1)
    out <- modelGeneVar(dummy, design=design0)
    ref <- modelGeneVar(dummy)
    expect_equal(out, ref)

    REFFUN <- function(y, X) {
        refit <- lm.fit(y=t(y), x=X)
        effects <- refit$effects[-seq_len(ncol(X)),]
        colMeans(effects^2)
    }

    # Testing with a one-way design matrix.
    design <- model.matrix(~factor(rep(c(1,2), each=100)))
    out2 <- modelGeneVar(dummy, design=design)
    expect_equal(ref$mean, out2$mean)

    test.var <- REFFUN(dummy, design)
    expect_equivalent(out2$total, test.var)

    # Testing with a more complex design matrix.
    Y <- runif(ncol(dummy))
    design <- model.matrix(~Y)
    out3 <- modelGeneVar(dummy, design=design)
    test.var <- REFFUN(dummy, design)

    expect_equal(ref$mean, out3$mean)
    expect_equivalent(out3$total, test.var)
})

test_that("modelGeneVar works with design and blocking", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    Y <- runif(ncol(dummy))
    design <- model.matrix(~Y)
    out <- modelGeneVar(dummy, design, block=block)

    accumulated.mean <- accumulated.total <- accumulated.tech <- 0
    for (i in unique(block)) {
        chosen <- i==block
        ref <- modelGeneVar(dummy[,chosen], design=design[chosen,])
        subout <- out$per.block[[i]]

        expect_identical(ref$mean, subout$mean)
        expect_identical(ref$total, subout$total)
        expect_identical(ref$tech, subout$tech)
        expect_identical(ref$bio, subout$bio)
        expect_identical(ref$p.value, subout$p.value)
    }
})

test_that("modelGeneVar works with SingleCellExperiment objects", {
    X <- SingleCellExperiment(list(logcounts=dummy))
    expect_equal(modelGeneVar(X), modelGeneVar(dummy))

    X <- SingleCellExperiment(list(whee=dummy))
    expect_equal(modelGeneVar(X, assay.type="whee"), modelGeneVar(dummy))
})
