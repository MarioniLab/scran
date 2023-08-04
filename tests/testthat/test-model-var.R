# This tests the functionality of modelGeneVar.
# library(testthat); library(scran); source("test-model-var.R")

set.seed(20001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

library(scuttle)
dummy <- normalizeCounts(dummy)

test_that("modelGeneVar works correctly without blocking", {
    out <- modelGeneVar(dummy)
    expect_equal(out$mean, unname(rowMeans(dummy)))
    expect_equal(out$total, unname(DelayedMatrixStats::rowVars(dummy)))
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
    expect_equal(out$p.value, metapod::parallelFisher(all.p)$p.value)

    # Responds to choice of method. 
    out2 <- modelGeneVar(dummy, block=block, method="stouffer")
    all.p <- lapply(out2$per.block, "[[", i="p.value")
    expect_equal(out2$p.value, metapod::parallelStouffer(all.p)$p.value)
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
    expect_equal(out$p.value, metapod::parallelFisher(all.p)$p.value)

    # Responds to choice of method with weighting.
    out2 <- modelGeneVar(dummy, block=block, method="stouffer", equiweight=FALSE)
    all.p <- lapply(out2$per.block, "[[", i="p.value")
    w <- countMatches(names(all.p), block)
    expect_equal(out2$p.value, metapod::parallelStouffer(all.p, weights=w)$p.value)
})

test_that("modelGeneVar handles blocks with no residual d.f.", {
    out <- modelGeneVar(dummy, block=rep(1:2, c(1, ncells-1)))
    ref <- modelGeneVar(dummy[,-1])
    expect_identical(out$mean, ref$mean)
    expect_identical(out$total, ref$total)

    out <- modelGeneVar(dummy, block=rep(1:3, c(1, 1, ncells-2)))
    ref <- modelGeneVar(dummy[,-c(1,2)])
    expect_identical(out$mean, ref$mean)
    expect_identical(out$total, ref$total)

    expect_error(modelGeneVar(dummy[,1,drop=FALSE]), "no residual d.f. in any level")

    # Also handles blocks with no observations whatsoever.
    block <- sample(2:5, ncells, replace=TRUE)
    ref  <- modelGeneVar(dummy, block=block)
    ref.empty <- modelGeneVar(dummy, block=factor(block, seq_len(5)))

    empty <- ref.empty$per.block[,"1"]
    expect_true(all(is.na(empty$p.value)))
    ref.empty$per.block <- ref.empty$per.block[,-1]
    expect_equal(ref.empty, ref)
})

test_that("modelGeneVar handles NAs in the blocking factor", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    keep <- !block %in% "A"
    ref <- modelGeneVar(dummy[,keep], block=block[keep])
    
    block2 <- block
    block2[!keep] <- NA
    out <- modelGeneVar(dummy, block=block2)

    expect_equal(out, ref)
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
    metadata(out) <- metadata(ref) <- list() # the fit parameters are more sensitive to numerical precision.
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

    # Flips out correctly.
    expect_error(modelGeneVar(dummy, design, block=sample(LETTERS[1:3], ncol(dummy), replace=TRUE)), "cannot specify 'design'")
    expect_error(modelGeneVar(dummy, design=diag(ncol(design))), "no residual d.f.")
})

test_that("modelGeneVar works with SingleCellExperiment objects", {
    X <- SingleCellExperiment(list(logcounts=dummy))
    expect_equal(modelGeneVar(X), modelGeneVar(dummy))

    X <- SingleCellExperiment(list(whee=dummy))
    expect_equal(modelGeneVar(X, assay.type="whee"), modelGeneVar(dummy))
})

test_that("modelGeneVar works with sparse inputs", {
    X <- dummy
    X_ <- as(X, "dgCMatrix")
    expect_equal(modelGeneVar(X), modelGeneVar(X_))

    # Safe with ultra-sparse rows.
    X[1:20,] <- 0
    X[1,2] <- 1
    X[2,c(10, 20, 30)] <- 1
    X[3,c(10, 20, 30)] <- 3:1
    X_ <- as(X, "dgCMatrix")
    expect_equal(modelGeneVar(X), modelGeneVar(X_))

    # Works with ultra dense rows.
    X <- dummy + 20
    X_ <- as(X, "dgCMatrix")
    expect_equal(modelGeneVar(X), modelGeneVar(X_))
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
    ref <- modelGeneVar(scuttle::normalizeCounts(dummy))
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    lspikes <- scuttle::normalizeCounts(spikes)
    expect_equal(metadata(out)$mean, rowMeans(lspikes))
    expect_equal(unname(metadata(out)$var), unname(DelayedMatrixStats::rowVars(lspikes)))

    fit <- fitTrendVar(metadata(out)$mean, metadata(out)$var)
    expect_identical(fit$std.dev, metadata(out)$std.dev)

    expect_equal(out$tech, fit$trend(ref$mean))
    expect_equal(out$bio, out$total - out$tech)
})

test_that("modelGeneVarWithSpikes works correctly with blocking", {
    block <- sample(LETTERS[1:5], ncells, replace=TRUE)
    out <- modelGeneVarWithSpikes(dummy, spikes, block=block)

    ref <- modelGeneVar(scuttle::normalizeCounts(dummy), block=block)
    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    accumulated.mean <- accumulated.total <- accumulated.tech <- 0
    sf1 <- scuttle::librarySizeFactors(dummy)
    sf2 <- scuttle::librarySizeFactors(spikes)

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

        # Increase 'tol' to avoid weird errors on Windows 32, 
        # probably due to FMA operations on the 32-bit x86 CRT. 
        expect_equal(ref$mean, subout$mean, tol=1e-6)
        expect_equal(ref$total, subout$total, tol=1e-6)
        expect_equal(ref$tech, subout$tech, tol=1e-6)
        expect_equal(ref$bio, subout$bio, tol=1e-6)
        expect_equal(ref$p.value, subout$p.value, tol=1e-6)

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
    expect_equal(out$p.value, metapod::parallelFisher(all.p)$p.value)
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
        ssf2 <- sf2[current] # use of 'sf2' is deliberate, avoid test errors due to slight numerical imprecision.
        ssf2 <- ssf2/mean(ssf2) * mean(ssf1)

        ref <- modelGeneVarWithSpikes(t(t(dummy[,current])/ssf1),
            size.factors=rep(1, sum(current)), 
            spikes=t(t(spikes[,current])/ssf2),
            spike.size.factors=rep(1, sum(current)))
        subout <- out$per.block[[i]]

        expect_equal(ref$mean, subout$mean)
        expect_equal(ref$total, subout$total)

        # Weird error propagation
        expect_equal(ref$tech, subout$tech, tol=1e-7)
        expect_equal(ref$bio, subout$bio, tol=1e-7)
        expect_equal(ref$p.value, subout$p.value, tol=1e-7)
    }
})

test_that("modelGeneVarWithSpikes works with design matrices", {
    Y <- runif(ncol(dummy))
    design <- model.matrix(~Y)

    genes <- modelGeneVar(scuttle::normalizeCounts(dummy), design=design)
    spiked <- modelGeneVar(scuttle::normalizeCounts(spikes), design=design)
    out <- modelGeneVarWithSpikes(dummy, spikes, design=design)

    expect_equal(out$mean, genes$mean)
    expect_equal(out$total, genes$total)
    expect_equal(metadata(out)$mean, setNames(spiked$mean, rownames(spikes)))
    expect_equal(metadata(out)$var, setNames(spiked$total, rownames(spikes)))
})

test_that("modelGeneVarWith Spikesworks with SingleCellExperiment objects", {
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
    expect_equal(modelGeneVarWithSpikes(X, "spikes", assay.type="whee"), 
        modelGeneVarWithSpikes(dummy, size.factors=sf1, spikes, spike.size.factors=sf2))
})

test_that("modelGeneVarWithSpikes works with a different pseudo-count", {
    out <- modelGeneVarWithSpikes(dummy, spikes, pseudo.count=3)
    ref <- modelGeneVar(scuttle::normalizeCounts(dummy, pseudo.count=3))

    expect_equal(out$mean, ref$mean)
    expect_equal(out$total, ref$total)

    lspikes <- scuttle::normalizeCounts(spikes, pseudo.count=3)
    expect_equal(metadata(out)$mean, rowMeans(lspikes))
    expect_equal(unname(metadata(out)$var), unname(DelayedMatrixStats::rowVars(lspikes)))
})

test_that("modelGeneVarWithSpikes works with sparse inputs", {
    ref <- modelGeneVarWithSpikes(dummy, spikes)
    d2 <- as(dummy, "dgCMatrix")
    s2 <- as(spikes, "dgCMatrix")
    expect_equal(ref, modelGeneVarWithSpikes(dummy, s2))
    expect_equal(ref, modelGeneVarWithSpikes(d2, spikes))
    expect_equal(ref, modelGeneVarWithSpikes(d2, s2))

    ref2 <- modelGeneVarWithSpikes(dummy, spikes, pseudo.count=0.5)
    expect_equal(ref2, modelGeneVarWithSpikes(d2, s2, pseudo.count=0.5))
})
