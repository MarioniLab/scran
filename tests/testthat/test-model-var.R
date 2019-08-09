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

    # Flips out if the design matrix is specified with block.
    expect_error(modelGeneVar(dummy, design, block=sample(LETTERS[1:3], ncol(dummy), replace=TRUE)), "cannot specify 'design'")
})

test_that("modelGeneVar works with SingleCellExperiment objects", {
    X <- SingleCellExperiment(list(logcounts=dummy))
    expect_equal(modelGeneVar(X), modelGeneVar(dummy))

    X <- SingleCellExperiment(list(whee=dummy))
    expect_equal(modelGeneVar(X, assay.type="whee"), modelGeneVar(dummy))
})

#######################################
#######################################
#######################################

set.seed(201001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

nspikes <- 100
smeans <- 2^runif(nspikes, -1, 5)
spikes <- matrix(rnbinom(nspikes*ncells, mu=smeans, size=5), ncol=ncells, nrow=nspikes)
rownames(spikes) <- paste0("X", seq_len(nspikes))

test_that("modelGeneVarWithSpikes works correctly in the basic case", {
    out <- modelGeneVarWithSpikes(dummy, spikes)
    ref <- modelGeneVar(scater::normalizeCounts(dummy))
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    lspikes <- scater::normalizeCounts(spikes)
    expect_equal(metadata(out)$mean, rowMeans(lspikes))
    expect_equal(unname(metadata(out)$var), DelayedMatrixStats::rowVars(lspikes))

    fit <- fitTrendVar(metadata(out)$mean, metadata(out)$var)
    expect_identical(fit$std.dev, metadata(out)$std.dev)

    expect_equal(out$tech, fit$trend(ref$mean))
    expect_equal(out$bio, out$total - out$tech)
})

test_that("modelGeneVarWithSpikes works correctly with blocking", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneVarWithSpikes(dummy, spikes, block=block)

    ref <- modelGeneVar(scater::normalizeCounts(dummy), block=block)
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    accumulated.mean <- accumulated.total <- accumulated.tech <- 0
    sf1 <- scater::librarySizeFactors(dummy)
    sf2 <- scater::librarySizeFactors(spikes)

    for (i in unique(block)) {
        current <- i==block

        # Forcibly avoid auto-centering of size.factors, to use as a reference here. 
        ssf1 <- sf1[current]
        ssf2 <- sf2[current]
        ssf2 <- ssf2/mean(ssf2) * mean(ssf1)

        ref <- modelGeneVarWithSpikes(t(t(dummy[,current])/ssf1),
            size.factors=rep(1, sum(current)), 
            spikes=t(t(spikes[,current])/ssf2),
            spike.size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_equal(ref$mean, subout$mean)
        expect_equal(ref$total, subout$total)
        expect_equal(ref$tech, subout$tech)
        expect_equal(ref$bio, subout$bio)
        expect_equal(ref$p.value, subout$p.value)

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
})

test_that("modelGeneVarWithSpikes centers size factors correctly", {
    # Without blocking.
    sf1 <- 2^rnorm(ncells, 0.05)
    sf2 <- 2^rnorm(ncells, 0.05)
    out <- modelGeneVarWithSpikes(dummy, size.factors=sf1, spikes=spikes, spike.size.factors=sf2)

    msf1 <- sf1/mean(sf1)
    msf2 <- sf2/mean(sf2)
    ref <- modelGeneVarWithSpikes(t(t(dummy)/msf1), size.factors=rep(1, ncells), 
        spikes=t(t(spikes)/msf2), spike.size.factors=rep(1, ncells))

    expect_equal(ref$mean, out$mean)
    expect_equal(ref$total, out$total)
    expect_equal(ref$tech, out$tech)
    expect_equal(ref$bio, out$bio)
    expect_equal(ref$p.value, out$p.value)

    # With blocking.
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneVarWithSpikes(dummy, size.factors=sf1, spikes=spikes, spike.size.factors=sf2, block=block)

    for (i in unique(block)) {
        current <- i==block

        ssf1 <- msf1[current]
        ssf2 <- msf2[current]
        ssf2 <- ssf2/mean(ssf2) * mean(ssf1)

        ref <- modelGeneVarWithSpikes(t(t(dummy[,current])/ssf1),
            size.factors=rep(1, sum(current)), 
            spikes=t(t(spikes[,current])/ssf2),
            spike.size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_equal(ref$mean, subout$mean)
        expect_equal(ref$total, subout$total)
        expect_equal(ref$tech, subout$tech)
        expect_equal(ref$bio, subout$bio)
        expect_equal(ref$p.value, subout$p.value)
    }
})

test_that("modelGeneVarWithSpikes works with design matrices", {
    Y <- runif(ncol(dummy))
    design <- model.matrix(~Y)

    genes <- modelGeneVar(scater::normalizeCounts(dummy), design=design)
    spiked <- modelGeneVar(scater::normalizeCounts(spikes), design=design)
    out <- modelGeneVarWithSpikes(dummy, spikes, design=design)

    expect_equal(out$mean, genes$mean)
    expect_equal(out$total, genes$total)
    expect_equal(metadata(out)$mean, spiked$mean)
    expect_equal(metadata(out)$var, spiked$total)
})

test_that("modelGeneVar works with SingleCellExperiment objects", {
    X <- SingleCellExperiment(list(counts=dummy))
    altExp(X, "spikes") <- SingleCellExperiment(list(counts=spikes))
    expect_equal(modelGeneVarWithSpikes(X, spikes="spikes"), modelGeneVarWithSpikes(dummy, spikes))

    X <- SingleCellExperiment(list(whee=dummy))
    altExp(X, "spikes") <- SingleCellExperiment(list(whee=spikes))
    expect_equal(modelGeneVarWithSpikes(X, "spikes", assay.type="whee"), modelGeneVarWithSpikes(dummy, spikes))

    X <- SingleCellExperiment(list(whee=dummy))
    sizeFactors(X) <- sf1 <- 2^rnorm(ncells, 0.1)
    altExp(X, "spikes") <- SingleCellExperiment(list(whee=spikes))
    sizeFactors(altExp(X)) <- sf2 <- 2^rnorm(ncells, 0.1)
    expect_equal(modelGeneVarWithSpikes(X, "spikes", assay.type="whee"), modelGeneVarWithSpikes(dummy, size.factors=sf1, spikes, spike.size.factors=sf2))
})
