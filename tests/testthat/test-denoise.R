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
    # Checking that the filtering for positive bio.comp and calculation of the variance explained is correct.
    npcs <- denoisePCA(lcounts, technical=fit$trend, value="n")
    keep <- rownames(dec)[dec$bio > 0]
    pc.out <- prcomp(t(lcounts[keep,]))

    tech.var <- fit$trend(rowMeans(lcounts))
    total.tech <- sum(tech.var[keep])
    expect_equal(npcs, scran:::.get_npcs_to_keep(pc.out$sdev^2, total.tech))

    expect_true(sum(pc.out$sdev[(npcs+1):ncol(lcounts)]^2) + pc.out$sdev[npcs+1]^2 * npcs < total.tech)
    expect_true(sum(pc.out$sdev[(npcs):ncol(lcounts)]^2) + pc.out$sdev[npcs]^2 * (npcs - 1) > total.tech)

    # Checking with different values for the technical noise, just in case.
    lower.fun <- function(x) { fit$trend(x) - 0.1 }
    npcs2 <- denoisePCA(lcounts, technical=lower.fun, value="n")
    expect_false(npcs==npcs2)

    tech.var2 <- lower.fun(rowMeans(lcounts))
    keep2 <- apply(lcounts, 1, var) > tech.var2
    pc.out2 <- prcomp(t(lcounts[keep2,]))
    total.tech2 <- sum(tech.var2[keep2])
    expect_equal(npcs2, scran:::.get_npcs_to_keep(pc.out2$sdev^2, total.tech2))

    expect_true(sum(pc.out2$sdev[(npcs2+1):ncol(lcounts)]^2) + pc.out2$sdev[npcs2+1]^2 * npcs2 < total.tech2)
    expect_true(sum(pc.out2$sdev[(npcs2):ncol(lcounts)]^2) + pc.out2$sdev[npcs2]^2 * (npcs2 - 1) > total.tech2)

    # And again.
    even.lower.fun <- function(x) { fit$trend(x) - 0.2 }
    npcs3 <- denoisePCA(lcounts, technical=even.lower.fun, value="n")
    expect_false(npcs==npcs3)
    
    tech.var3 <- even.lower.fun(rowMeans(lcounts))
    keep3 <- apply(lcounts, 1, var) > tech.var3
    pc.out3 <- prcomp(t(lcounts[keep3,]))
    total.tech3 <- sum(tech.var3[keep3])
    expect_equal(npcs3, scran:::.get_npcs_to_keep(pc.out3$sdev^2, total.tech3))

    expect_true(sum(pc.out3$sdev[(npcs3+1):ncol(lcounts)]^2) + pc.out3$sdev[npcs3+1]^2 * npcs3 < total.tech3)
    expect_true(sum(pc.out3$sdev[(npcs3):ncol(lcounts)]^2) + pc.out3$sdev[npcs3]^2 * (npcs3 - 1) > total.tech3)

    # Checking that the PC selection is correct.
    pcs <- denoisePCA(lcounts, technical=fit$trend, value="pca")
    expect_equal(pcs, pc.out$x[,seq_len(npcs)])
    pcs2 <- denoisePCA(lcounts, technical=lower.fun, value="pca")
    expect_equal(pcs2, pc.out2$x[,seq_len(npcs2)])
    pcs3 <- denoisePCA(lcounts, technical=even.lower.fun, value="pca")
    expect_equal(pcs3, pc.out3$x[,seq_len(npcs3)])

    # Checking that the low-rank approximation is correctly computed.
    lrout <- denoisePCA(lcounts, technical=fit$trend, value="lowrank")
    expect_identical(dim(lrout), dim(lcounts))
    expect_equal(rowMeans(lrout)[keep], rowMeans(lcounts)[keep])

    expect_true(all(lrout[setdiff(rownames(lrout), keep),]==0)) # internally filtered genes set to all-zero.
    expect_true(all(apply(lrout[keep,], 1, var) > 0))

    QR <- qr(lrout - rowMeans(lrout)) 
    expect_equal(QR$rank, npcs) # checking that it has the correct rank.
    expect_equal(sum(apply(pcs, 2, var)), sum(apply(lrout, 1, var))) # explains the same amount of variance.
}) 

test_that("denoisePCA works with different settings", {
    # Checking that the output is the same.
    pcs <- denoisePCA(lcounts, technical=fit$trend)
    pcs3 <- denoisePCA(lcounts, technical=fit$trend, design=cbind(rep(1, ncells)))
    are_PCs_equal(pcs, pcs3)

    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
    pcs2 <- denoisePCA(lcounts[not.spike,], technical=fit$trend)
    are_PCs_equal(pcs, pcs2)

    # Checking that low rank settings behave correctly with subsetting.
    lr1 <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike, value="lowrank")
    lr2 <- denoisePCA(lcounts[not.spike,], technical=fit$trend, value="lowrank")
    expect_equal(lr1, lr2)

    lr3 <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike, value="lowrank", preserve.dim=TRUE)
    expect_equal(lr1, lr3[not.spike,])
    expect_true(all(lr3[-not.spike,]==0))

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
    
    keep <- true.var > dfit$trend(rowMeans(lcounts))
    alt.pc <- prcomp(new.x[,keep])
    are_PCs_equal(alt.pc$x[,seq_len(ncol(pcs)),drop=FALSE], pcs)
})

test_that("denoisePCA throws errors correctly", {
    # Checking invalid specifications.
    expect_error(denoisePCA(lcounts[0,], fit$trend), "a dimension is zero")
    expect_error(denoisePCA(lcounts[,0], fit$trend), "error code")
})

test_that("denoisePCA works with SingleCellExperiment inputs", {
    # Checking for proper behaviour with SCESet.
    X <- SingleCellExperiment(list(exprs=lcounts))
    X2 <- denoisePCA(X, technical=fit$trend)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    pcs <- denoisePCA(lcounts, technical=fit$trend)
    are_PCs_equal(pcx, pcs)
    
    isSpike(X, "Spike") <- is.spike
    X2 <- denoisePCA(X, technical=fit$trend, get.spikes=TRUE)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    are_PCs_equal(pcx, pcs)
    
    X2 <- denoisePCA(X, technical=fit$trend)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
    are_PCs_equal(pcx, pcs)

    # Checking lowrank calculations.
    X3 <- denoisePCA(X, technical=fit$trend, value="lowrank")
    ref <- denoisePCA(exprs(X), technical=fit$trend, value="lowrank", subset.row=not.spike)
    pcx <- assay(X3, "lowrank")
    expect_equal(pcx[not.spike,], ref)
    expect_true(all(pcx[is.spike,]==0))
    
    X3 <- denoisePCA(X, technical=fit$trend, value="lowrank", get.spikes=TRUE)
    ref <- denoisePCA(exprs(X), technical=fit$trend, value="lowrank")
    pcx <- assay(X3, "lowrank")
    expect_equal(pcx, ref)
})
