# This tests the functionality of modelGeneCV2.
# library(testthat); library(scran); source("test-model-cv2.R")

set.seed(20001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

library(scuttle)
dummy2 <- normalizeCounts(dummy, log=FALSE)

test_that("modelGeneCV2 works correctly without blocking", {
    out <- modelGeneCV2(dummy)
    expect_equal(out$mean, unname(rowMeans(dummy2)))
    expect_equal(out$total, unname(DelayedMatrixStats::rowVars(dummy2)/out$mean^2))
    expect_equal(out$trend, metadata(out)$trend(out$mean))
    expect_equal(out$ratio, out$total/out$trend)
    expect_equal(order(out$p.value), order(out$ratio, decreasing=TRUE))
})

test_that("modelGeneCV2 responds to size factors", {
    sf <- runif(ncells)
    out <- modelGeneCV2(dummy, size.factors=sf)
    ref <- modelGeneCV2(t(t(dummy)/sf*mean(sf)), size.factors=rep(1, ncells))
    expect_equal(out, ref, tol = 1e-6)
})

test_that("modelGeneCV2 works correctly with blocking, no weighting", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneCV2(dummy, block=block)

    accumulated.mean <- accumulated.total <- accumulated.trend <- 0
    for (i in unique(block)) {
        current <- i==block
        ref <- modelGeneCV2(dummy2[,current], size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_identical(ref$mean, subout$mean)
        expect_identical(ref$total, subout$total)
        expect_identical(ref$trend, subout$trend)
        expect_identical(ref$ratio, subout$ratio)
        expect_identical(ref$p.value, subout$p.value)

        accumulated.mean <- accumulated.mean + log(ref$mean)
        accumulated.total <- accumulated.total + log(ref$total)
        accumulated.trend <- accumulated.trend + log(ref$trend)
    }

    # Check combining statistics works correctly.
    n <- length(unique(block))
    expect_equal(out$mean, exp(accumulated.mean/n))
    expect_equal(out$total, exp(accumulated.total/n))
    expect_equal(out$trend, exp(accumulated.trend/n))
    expect_equal(out$ratio, out$total/out$trend)

    all.p <- lapply(out$per.block, "[[", i="p.value")
    expect_equal(out$p.value, metapod::parallelFisher(all.p)$p.value)

    # Responds to choice of method. 
    out2 <- modelGeneCV2(dummy, block=block, method="stouffer")
    all.p <- lapply(out2$per.block, "[[", i="p.value")
    expect_equal(out2$p.value, metapod::parallelStouffer(all.p)$p.value)
})

test_that("modelGeneCV2 works correctly with blocking and weighting", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneCV2(dummy, block=block, equiweight=FALSE)

    accumulated.mean <- accumulated.total <- accumulated.trend <- 0
    for (i in unique(block)) {
        current <- i==block
        ref <- modelGeneCV2(dummy2[,current], size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_identical(ref$mean, subout$mean)
        expect_identical(ref$total, subout$total)
        expect_identical(ref$trend, subout$trend)
        expect_identical(ref$ratio, subout$ratio)
        expect_identical(ref$p.value, subout$p.value)

        n <- sum(i==block)
        accumulated.mean <- accumulated.mean + log(ref$mean) * n
        accumulated.total <- accumulated.total + log(ref$total) * n
        accumulated.trend <- accumulated.trend + log(ref$trend) * n
    }

    # Check combining statistics works correctly.
    n <- length(block)
    expect_equal(out$mean, exp(accumulated.mean/n))
    expect_equal(out$total, exp(accumulated.total/n))
    expect_equal(out$trend, exp(accumulated.trend/n))
    expect_equal(out$ratio, out$total/out$trend)

    all.p <- lapply(out$per.block, "[[", i="p.value")
    expect_equal(out$p.value, metapod::parallelFisher(all.p)$p.value)

    # Responds to choice of method with weighting.
    out2 <- modelGeneCV2(dummy, block=block, method="stouffer", equiweight=FALSE)
    all.p <- lapply(out2$per.block, "[[", i="p.value")
    w <- countMatches(names(all.p), block)
    expect_equal(out2$p.value, metapod::parallelStouffer(all.p, weights=w)$p.value)
})

test_that("modelGeneCV2 handles blocks with no residual d.f.", {
    out <- modelGeneCV2(dummy2, size.factors=rep(1, ncells), block=rep(1:2, c(1, ncells-1)))
    ref <- modelGeneCV2(dummy2[,-1], size.factors=rep(1, ncells))
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    out <- modelGeneCV2(dummy, size.factors=rep(1, ncells), block=rep(1:3, c(1, 1, ncells-2)))
    ref <- modelGeneCV2(dummy[,-c(1,2)], size.factors=rep(1, ncells))
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    expect_error(modelGeneCV2(dummy[,1,drop=FALSE]), "no residual d.f. in any level")
})

test_that("modelGeneCV2 works with subsetting options", {
    chosen <- sample(ngenes, ngenes/2)
    out <- modelGeneCV2(dummy, subset.row=chosen)
    ref <- modelGeneCV2(dummy[chosen,])
    expect_equal(out, ref)

    # Subsetting by fit works.
    out2 <- modelGeneCV2(dummy, subset.fit=chosen, size.factors=librarySizeFactors(dummy, chosen))
    expect_identical(rownames(out2), rownames(dummy))
    expect_equal(out2[chosen,1:5], ref[,1:5])

    # Zero-length subsetting.
    empty <- modelGeneCV2(dummy, subset.row=integer(0), subset.fit=chosen, size.factors=librarySizeFactors(dummy))
    expect_identical(nrow(empty), 0L)

    expect_error(modelGeneCV2(dummy, subset.fit=integer(0)), "need at least 2 points")
})

test_that("modelGeneCV2 works with SingleCellExperiment objects", {
    X <- SingleCellExperiment(list(counts=dummy))
    expect_equal(modelGeneCV2(X), modelGeneCV2(dummy))

    sizeFactors(X) <- runif(ncol(X), 0.5, 1.5)
    expect_equal(modelGeneCV2(X), modelGeneCV2(dummy, sizeFactors(X)))

    X <- SingleCellExperiment(list(whee=dummy))
    expect_equal(modelGeneCV2(X, assay.type="whee"), modelGeneCV2(dummy))
})

test_that("modelGeneCV2 works with sparse inputs", {
    X <- dummy2
    sf <- librarySizeFactors(X) # need to specify this for testing, as minor differences in colSums propagate.
    X_ <- as(X, "dgCMatrix")
    expect_equal(modelGeneCV2(X, size.factors=sf), modelGeneCV2(X_, size.factors=sf))

    # Safe with ultra-sparse rows.
    X[1:20,] <- 0
    X[1,2] <- 1
    X[2,c(10, 20, 30)] <- 1
    X[3,c(10, 20, 30)] <- 3:1
    X_ <- as(X, "dgCMatrix")
    expect_equal(modelGeneCV2(X, size.factors=sf), modelGeneCV2(X_, size.factors=sf))

    # Works with ultra dense rows.
    X <- dummy2 + 1
    X_ <- as(X, "dgCMatrix")
    expect_equal(modelGeneCV2(X, size.factors=sf), modelGeneCV2(X_, size.factors=sf))
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

normdummy <- scuttle::normalizeCounts(dummy, log=FALSE)
normspikes <- scuttle::normalizeCounts(spikes, log=FALSE)

test_that("modelGeneCV2WithSpikes works correctly in the basic case", {
    out <- modelGeneCV2WithSpikes(dummy, spikes)
    ref <- modelGeneCV2(normdummy)
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    expect_equal(metadata(out)$mean, rowMeans(normspikes))
    expect_equal(metadata(out)$cv2, DelayedMatrixStats::rowVars(normspikes)/metadata(out)$mean^2)

    fit <- fitTrendCV2(metadata(out)$mean, metadata(out)$cv2, ncells)
    expect_identical(fit$std.dev, metadata(out)$std.dev)

    expect_equal(out$trend, fit$trend(ref$mean))
    expect_equal(out$ratio, out$total/out$trend)
})

test_that("modelGeneCV2WithSpikes works correctly with blocking", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneCV2WithSpikes(dummy, spikes, block=block)

    ref <- modelGeneCV2(normdummy, block=block)
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    accumulated.mean <- accumulated.total <- accumulated.trend <- 0
    sf1 <- scuttle::librarySizeFactors(dummy)
    sf2 <- scuttle::librarySizeFactors(spikes)

    for (i in unique(block)) {
        current <- i==block

        # Forcibly avoid auto-centering of size.factors, to use as a reference here. 
        ssf1 <- sf1[current]
        ssf2 <- sf2[current]
        ssf2 <- ssf2/mean(ssf2) * mean(ssf1)

        ref <- modelGeneCV2WithSpikes(t(t(dummy[,current])/ssf1),
            size.factors=rep(1, sum(current)), 
            spikes=t(t(spikes[,current])/ssf2),
            spike.size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_equal(ref$mean, subout$mean)
        expect_equal(ref$total, subout$total)
        expect_equal(ref$trend, subout$trend)
        expect_equal(ref$ratio, subout$ratio)
        expect_equal(ref$p.value, subout$p.value)

        accumulated.mean <- accumulated.mean + log(ref$mean)
        accumulated.total <- accumulated.total + log(ref$total)
        accumulated.trend <- accumulated.trend + log(ref$trend)
    }

    # Check combining statistics works correctly.
    n <- length(unique(block))
    expect_equal(out$mean, exp(accumulated.mean/n))
    expect_equal(out$total, exp(accumulated.total/n))
    expect_equal(out$trend, exp(accumulated.trend/n))
    expect_equal(out$ratio, out$total/out$trend)

    all.p <- lapply(out$per.block, "[[", i="p.value")
    expect_equal(out$p.value, metapod::parallelFisher(all.p)$p.value)
})

test_that("modelGeneCV2WithSpikes centers size factors correctly", {
    # Without blocking.
    sf1 <- 2^rnorm(ncells, 0.05)
    sf2 <- 2^rnorm(ncells, 0.05)
    out <- modelGeneCV2WithSpikes(dummy, size.factors=sf1, spikes=spikes, spike.size.factors=sf2)

    msf1 <- sf1/mean(sf1)
    msf2 <- sf2/mean(sf2)
    ref <- modelGeneCV2WithSpikes(t(t(dummy)/msf1), size.factors=rep(1, ncells), 
        spikes=t(t(spikes)/msf2), spike.size.factors=rep(1, ncells))

    expect_equal(ref$mean, out$mean, tol=1e-6)
    expect_equal(ref$total, out$total, tol=1e-6)
    expect_equal(ref$trend, out$trend, tol=1e-6)
    expect_equal(ref$ratio, out$ratio, tol=1e-6)
    expect_equal(ref$p.value, out$p.value, tol=1e-6)

    # With blocking.
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneCV2WithSpikes(dummy, size.factors=sf1, spikes=spikes, spike.size.factors=sf2, block=block)

    for (i in unique(block)) {
        current <- i==block

        ssf1 <- msf1[current]
        ssf2 <- msf2[current]
        ssf2 <- ssf2/mean(ssf2) * mean(ssf1)

        ref <- modelGeneCV2WithSpikes(t(t(dummy[,current])/ssf1),
            size.factors=rep(1, sum(current)), 
            spikes=t(t(spikes[,current])/ssf2),
            spike.size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_equal(ref$mean, subout$mean, tol=1e-6)
        expect_equal(ref$total, subout$total, tol=1e-6)
        expect_equal(ref$trend, subout$trend, tol=1e-6)
        expect_equal(ref$ratio, subout$ratio, tol=1e-6)
        expect_equal(ref$p.value, subout$p.value, tol=1e-6)
    }
})

test_that("modelGeneCV2 works with SingleCellExperiment objects", {
    X <- SingleCellExperiment(list(counts=dummy))
    altExp(X, "spikes") <- SingleCellExperiment(list(counts=spikes))
    expect_equal(modelGeneCV2WithSpikes(X, spikes="spikes"), modelGeneCV2WithSpikes(dummy, spikes))

    X <- SingleCellExperiment(list(whee=dummy))
    altExp(X, "spikes") <- SingleCellExperiment(list(whee=spikes))
    expect_equal(modelGeneCV2WithSpikes(X, "spikes", assay.type="whee"), modelGeneCV2WithSpikes(dummy, spikes))

    X <- SingleCellExperiment(list(whee=dummy))
    sizeFactors(X) <- sf1 <- 2^rnorm(ncells, 0.1)
    altExp(X, "spikes") <- SingleCellExperiment(list(whee=spikes))
    sizeFactors(altExp(X)) <- sf2 <- 2^rnorm(ncells, 0.1)
    expect_equal(modelGeneCV2WithSpikes(X, "spikes", assay.type="whee"), 
        modelGeneCV2WithSpikes(dummy, size.factors=sf1, spikes, spike.size.factors=sf2))
})

test_that("modelGeneCV2WithSpikes works with sparse inputs", {
    # Need to specify the size factors to avoid propagating small errors due to colSums differences.
    sf1 <- librarySizeFactors(normdummy)
    sf2 <- librarySizeFactors(normspikes)
    ref <- modelGeneVarWithSpikes(normdummy, normspikes, size.factors=sf1, spike.size.factors=sf2)

    d2 <- as(normdummy, "dgCMatrix")
    s2 <- as(normspikes, "dgCMatrix")
    expect_equal(ref, modelGeneVarWithSpikes(normdummy, s2, size.factors=sf1, spike.size.factors=sf2))
    expect_equal(ref, modelGeneVarWithSpikes(d2, normspikes, size.factors=sf1, spike.size.factors=sf2))
    expect_equal(ref, modelGeneVarWithSpikes(d2, s2, size.factors=sf1, spike.size.factors=sf2))
})
