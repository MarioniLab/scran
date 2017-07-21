# This checks the denoisePCA function.
# require(scran); require(testthat); source("test-denoise.R")

are_PCs_equal <- function(first, second) {
    expect_identical(dim(first), dim(second))
    relative <- first/second
    diffs <- abs(colSums(relative))
    expect_true(all(abs(diffs - nrow(first)) < 1e-8))
}

# Mocking up some data with subpopulations of cells.

set.seed(1000)
ngenes <- 1000
npops <- 5
ncells <- 100
means <- 2^runif(ngenes, -1, 10)
pops <- matrix(2^rnorm(npops * ngenes), ncol=npops) * means

is.spike <- 1:100
pops[is.spike,] <- means[is.spike] # spike ins are constant across subpopulations.
in.pop <- sample(npops, ncells, replace=TRUE)
true.means <- pops[,in.pop,drop=FALSE]

dispersions <- 10/means + 0.2
ncells <- 100
counts <- matrix(rnbinom(ngenes*ncells, mu=true.means, size=1/dispersions), ncol=ncells)
rownames(counts) <- paste0("Gene", seq_len(ngenes))

lcounts <- log2(counts + 1)
fit <- trendVar(lcounts, subset.row=is.spike)
dec <- decomposeVar(lcounts, fit)

test_that("denoisePCA works as expected", {
    # Checking that the calculation of the variance explained is correct.
    total.tech <- fit$trend(rowMeans(lcounts))
    npcs <- denoisePCA(lcounts, technical=fit$trend, value="n")
    pc.out <- prcomp(t(lcounts))
    expect_equal(npcs, scran:::.get_npcs_to_keep(pc.out$sdev^2, sum(total.tech)))

    # Checking with different values for the technical noise, just in case.
    npcs2 <- denoisePCA(lcounts, technical=total.tech - 0.1, value="n")
    expect_false(npcs==npcs2)
    expect_equal(npcs2, scran:::.get_npcs_to_keep(pc.out$sdev^2, sum(total.tech - 0.1)))
    npcs3 <- denoisePCA(lcounts, technical=total.tech - 0.2, value="n")
    expect_false(npcs==npcs3)
    expect_equal(npcs3, scran:::.get_npcs_to_keep(pc.out$sdev^2, sum(total.tech - 0.2)))

    # Checking that the PC calculation is correct.
    pcs <- denoisePCA(lcounts, technical=fit$trend, value="pca")
    expect_equal(pcs, pc.out$x[,seq_len(npcs)])
    pcs2 <- denoisePCA(lcounts, technical=total.tech - 0.1, value="pca")
    expect_equal(pcs2, pc.out$x[,seq_len(npcs2)])
    pcs3 <- denoisePCA(lcounts, technical=total.tech - 0.2, value="pca")
    expect_equal(pcs3, pc.out$x[,seq_len(npcs3)])

    # Checking that the low-rank approximation is correctly computed.
    lrout <- denoisePCA(lcounts, technical=fit$trend, value="lowrank")
    expect_identical(dim(lrout), dim(lcounts))
    expect_equal(rowMeans(lrout), rowMeans(lcounts))
    QR <- qr(lrout - rowMeans(lcounts)) # has the correct rank.
    expect_equal(QR$rank, npcs)
    expect_equal(sum(apply(pcs, 2, var)), sum(apply(lrout, 1, var))) # explains the same amount of variance.
}) 

test_that("denoisePCA works with different settings", {
    # Checking that the output is the same.
    pcs <- denoisePCA(lcounts, technical=fit$trend)
    pcs2 <- denoisePCA(lcounts, technical=setNames(dec$tech, rownames(dec)))
    are_PCs_equal(pcs, pcs2)
   
    pcs3 <- denoisePCA(lcounts, technical=fit$trend, design=cbind(rep(1, ncells)))
    are_PCs_equal(pcs, pcs3)

    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
    pcs2 <- denoisePCA(lcounts[not.spike,], technical=fit$trend)
    are_PCs_equal(pcs, pcs2)

    # Checking that it responds correctly to min and max settings.
    ref <- denoisePCA(lcounts, technical=fit$trend)
    pcs <- denoisePCA(lcounts, technical=fit$trend, min.rank=ncol(ref)+1)
    expect_identical(ncol(pcs), ncol(ref)+1L)
    pcs <- denoisePCA(lcounts, technical=fit$trend, max.rank=ncol(ref)-1)
    expect_identical(pcs, ref[,-ncol(ref)])
})

test_that("denoisePCA works with design matrices", {
    # Checking for sensible handling of design matrices.
    design <- model.matrix(~factor(in.pop))
    dfit <- trendVar(lcounts, subset.row=is.spike, design=design)
    pcs <- denoisePCA(lcounts, design=design, technical=dfit$trend)
    
    alt <- lm.fit(y=t(lcounts), x=design)
    true.var <- colMeans(alt$effects[-seq_len(alt$rank),]^2)
    obs.var <- apply(alt$residuals, 2, var)
    new.x <- alt$residuals * sqrt(true.var/obs.var)
    alt.pc <- prcomp(new.x)
    are_PCs_equal(alt.pc$x[,seq_len(ncol(pcs)),drop=FALSE], pcs)
})

test_that("denoisePCA throws errors correctly", {
    # Checking invalid specifications.
    expect_error(denoisePCA(lcounts, technical=c(Whee=1)), "missing gene names in 'technical'")
    unnamed.lcounts <- lcounts
    rownames(unnamed.lcounts) <- NULL
    expect_error(denoisePCA(unnamed.lcounts, technical=c(Whee=1)), "rows of 'x' should be named with gene names")
})

test_that("denoisePCA works with SCESet inputs", {
    # Checking for proper behaviour with SCESet.
    X <- newSCESet(exprsData=lcounts, logExprsOffset=1, lowerDetectionLimit=0)
    X2 <- denoisePCA(X, technical=fit$trend)
    pcx <- reducedDimension(X2)
    rownames(pcx) <- NULL
    pcs <- denoisePCA(lcounts, technical=fit$trend)
    are_PCs_equal(pcx, pcs)
    
    X <- calculateQCMetrics(X, feature_controls=list(Spike=is.spike))
    setSpike(X) <- "Spike"
    X2 <- denoisePCA(X, technical=fit$trend, get.spikes=TRUE)
    pcx <- reducedDimension(X2)
    rownames(pcx) <- NULL
    are_PCs_equal(pcx, pcs)
    
    X2 <- denoisePCA(X, technical=fit$trend)
    pcx <- reducedDimension(X2)
    rownames(pcx) <- NULL
    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
    are_PCs_equal(pcx, pcs)

    # Checking lowrank calculations.
    X3 <- denoisePCA(X, technical=fit$trend, value="lowrank")
    ref <- denoisePCA(exprs(X), technical=fit$trend, value="lowrank", subset.row=not.spike)
    pcx <- assayDataElement(X3, "lowrank")
    expect_equal(pcx[not.spike,], ref)
    expect_true(all(is.na(pcx[is.spike,])))
    
    X3 <- denoisePCA(X, technical=fit$trend, value="lowrank", get.spikes=TRUE)
    ref <- denoisePCA(exprs(X), technical=fit$trend, value="lowrank")
    pcx <- assayDataElement(X3, "lowrank")
    expect_equal(pcx, ref)
})
