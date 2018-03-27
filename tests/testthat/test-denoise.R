# This checks the denoisePCA function.
# require(scran); require(testthat); source("test-denoise.R")

are_PCs_equal <- function(first, second, tol=1e-8) {
    expect_identical(dim(first), dim(second))
    relative <- first/second
    expect_true(all(colSums(relative > 0) %in% c(0, nrow(first))))
    expect_true(all(abs(abs(relative)-1) < tol))
}

# Mocking up some data with subpopulations of cells.

set.seed(1000)
ngenes <- 1000
npops <- 10
ncells <- 200
means <- 2^runif(ngenes, -1, 10)
pops <- matrix(2^rnorm(npops * ngenes), ncol=npops) * means

is.spike <- 1:100
pops[is.spike,] <- means[is.spike] # spike ins are constant across subpopulations.
in.pop <- sample(npops, ncells, replace=TRUE)
true.means <- pops[,in.pop,drop=FALSE]

dispersions <- 10/means + 0.2
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

    verify_npcs <- function(npcs, sdev, tech.total) {
        var.exp <- sdev^2
        expect_equal(npcs[1], scran:::.get_npcs_to_keep(var.exp, tech.total))
        
        # Chosen number of PCs should be at the technical threshold.
        expect_true(sum(var.exp[(npcs+1):ncol(lcounts)]) < tech.total) 
        expect_true(sum(var.exp[(npcs):ncol(lcounts)]) > tech.total)
    
        reported <- attr(npcs, "percentVar")
        exp.var <- sdev^2
        expect_equal(reported, exp.var[seq_along(reported)]/sum(exp.var))
    }

    tech.var <- fit$trend(rowMeans(lcounts))
    total.tech <- sum(tech.var[keep])
    verify_npcs(npcs, pc.out$sdev, total.tech)

    # Checking with different values for the technical noise, just in case.
    lower.fun <- function(x) { fit$trend(x) - 0.1 }
    npcs2 <- denoisePCA(lcounts, technical=lower.fun, value="n")
    expect_false(npcs==npcs2)

    tech.var2 <- lower.fun(rowMeans(lcounts))
    keep2 <- apply(lcounts, 1, var) > tech.var2
    pc.out2 <- prcomp(t(lcounts[keep2,]))
    total.tech2 <- sum(tech.var2[keep2])

    verify_npcs(npcs2, pc.out2$sdev, total.tech2)

    # And again.
    even.lower.fun <- function(x) { fit$trend(x) - 0.2 }
    npcs3 <- denoisePCA(lcounts, technical=even.lower.fun, value="n")
    expect_false(npcs==npcs3)
    
    tech.var3 <- even.lower.fun(rowMeans(lcounts))
    keep3 <- apply(lcounts, 1, var) > tech.var3
    pc.out3 <- prcomp(t(lcounts[keep3,]))
    total.tech3 <- sum(tech.var3[keep3])

    verify_npcs(npcs3, pc.out3$sdev, total.tech3)

    # Checking that the PC selection is correct.
    pcs <- denoisePCA(lcounts, technical=fit$trend, value="pca")
    expect_equal(pcs[,], pc.out$x[,seq_len(npcs)])
    expect_equal(attr(pcs, "percentVar"), attr(npcs, "percentVar"))

    pcs2 <- denoisePCA(lcounts, technical=lower.fun, value="pca")
    expect_equal(pcs2[,], pc.out2$x[,seq_len(npcs2)])
    expect_equal(attr(pcs2, "percentVar"), attr(npcs2, "percentVar"))

    pcs3 <- denoisePCA(lcounts, technical=even.lower.fun, value="pca")
    expect_equal(pcs3[,], pc.out3$x[,seq_len(npcs3)])
    expect_equal(attr(pcs3, "percentVar"), attr(npcs3, "percentVar"))
})

test_that("Low-rank approximations work as expected", {
    # Checking that the low-rank approximation is correctly computed.
    lrout <- denoisePCA(lcounts, technical=fit$trend, value="lowrank")
    expect_identical(dim(lrout), dim(lcounts))
    
    keep <- rownames(dec)[dec$bio > 0]
    expect_equal(rowMeans(lrout)[keep], rowMeans(lcounts)[keep])
    expect_true(all(apply(lrout[keep,], 1, var) > 0))

    QR <- qr(lrout - rowMeans(lrout)) 
    pcs <- denoisePCA(lcounts, technical=fit$trend, value="pca")
    expect_equal(ncol(pcs), QR$rank) # checking that it has the correct rank.
    expect_equal(sum(apply(pcs, 2, var)), sum(apply(lrout[keep,], 1, var))) # explains the same amount of variance.

    # Checking that projections work.
    lrout <- denoisePCA(lcounts, technical=fit$trend, value="lowrank")
    lcounts.extra <- rbind(lcounts, lcounts[1:10,])
    lrout.extra <- denoisePCA(lcounts.extra, technical=fit$trend, value="lowrank", subset.row=seq_len(nrow(lcounts)))
    expect_equal(lrout[,], lrout.extra[seq_len(nrow(lcounts)),])
    expect_equal(lrout[1:10,], lrout.extra[nrow(lcounts)+seq_len(10),])
}) 

test_that("denoisePCA works with different settings", {
    # Checking proper behaviour with subsetting.
    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
    pcs2 <- denoisePCA(lcounts[not.spike,], technical=fit$trend)
    are_PCs_equal(pcs, pcs2)

    lr1 <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike, value="lowrank")
    lr2 <- denoisePCA(lcounts[not.spike,], technical=fit$trend, value="lowrank")
    expect_equal(lr1[not.spike,], lr2[,])

    # Checking that it responds correctly to min and max settings.
    ref <- denoisePCA(lcounts, technical=fit$trend)
    pcs <- denoisePCA(lcounts, technical=fit$trend, min.rank=ncol(ref)+1)
    expect_identical(ncol(pcs), ncol(ref)+1L)
    expect_identical(pcs[,seq_len(ncol(ref))], ref[,]) 

    pcs <- denoisePCA(lcounts, technical=fit$trend, max.rank=ncol(ref)-1)
    expect_identical(ncol(pcs), ncol(ref)-1L)
    expect_identical(pcs[,], ref[,-ncol(ref)])
})

test_that("denoisePCA works with IRLBA", {
    # Checking choice of number of PCs.
    keep <- dec$bio > 0
    posbio <- lcounts[keep,]
    df0 <- ncol(posbio)-1
    max.cells <- formals(scran:::.denoisePCA)$max.rank # using the upper limit in denoisePCA.
    current <- t(posbio - rowMeans(posbio))
    
    npcs <- suppressWarnings(denoisePCA(lcounts, technical=fit$trend, value="n", approximate=TRUE, rand.seed=100))
    set.seed(100)
    e1 <- suppressWarnings(irlba::irlba(current, nu=0, nv=max.cells))

    total.var <- sum(dec$total[keep])
    expect_equal(npcs[1], scran:::.get_npcs_to_keep(e1$d^2/df0, sum(dec$tech[keep]), total=total.var))
    expect_equal(attr(npcs, "percentVar"), e1$d^2/df0/total.var)
    
    # Checking the actual PCs themselves.
    pca <- suppressWarnings(denoisePCA(lcounts, technical=fit$trend, value="pca", approximate=TRUE, rand.seed=200))
    set.seed(200)
    epc <- suppressWarnings(irlba::prcomp_irlba(current, n=max.cells, center=FALSE, scale.=FALSE))
    are_PCs_equal(pca, epc$x[,seq_len(npcs),drop=FALSE]) 
    expect_equal(attr(pca, "percentVar"), epc$sdev^2/total.var)
    
    # Checking the low-rank approximations.
    lr <- suppressWarnings(denoisePCA(lcounts, technical=fit$trend, value="lowrank", approximate=TRUE, rand.seed=300))

    set.seed(300)
    e2 <- suppressWarnings(irlba::irlba(current, nu=max.cells, nv=max.cells)) 
    lowrank <- e2$u[,1:npcs] %*% (e2$d[1:npcs] * t(e2$v[,1:npcs]))   

    unnamed.lr <- lr
    dimnames(unnamed.lr) <- NULL
    expect_equivalent(unnamed.lr[keep,], t(lowrank) + unname(rowMeans(posbio)))
})

test_that("denoisePCA throws errors correctly", {
    # Checking invalid specifications.
    expect_error(denoisePCA(lcounts[0,], fit$trend), "a dimension is zero")
    expect_error(denoisePCA(lcounts[,0], fit$trend), "a dimension is zero")
})

test_that("denoisePCA works with SingleCellExperiment inputs", {
    # Checking for proper behaviour with SCESet.
    X <- SingleCellExperiment(list(logcounts=lcounts))
    X2 <- denoisePCA(X, technical=fit$trend)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    pcs <- denoisePCA(lcounts, technical=fit$trend)
    are_PCs_equal(pcx, pcs)
    expect_identical(attr(pcx, "percentVar"), attr(pcs, "percentVar"))
    
    isSpike(X, "Spike") <- is.spike
    X2 <- denoisePCA(X, technical=fit$trend, get.spikes=TRUE)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    are_PCs_equal(pcx, pcs)
    expect_identical(attr(pcx, "percentVar"), attr(pcs, "percentVar"))
    
    X2 <- denoisePCA(X, technical=fit$trend)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
    are_PCs_equal(pcx, pcs)
    expect_identical(attr(pcx, "percentVar"), attr(pcs, "percentVar"))

    # Checking lowrank calculations.
    X3 <- denoisePCA(X, technical=fit$trend, value="lowrank")
    ref <- denoisePCA(exprs(X)[not.spike,], technical=fit$trend, value="lowrank")
    pcx <- assay(X3, "lowrank")
    expect_equal(pcx[not.spike,], ref[,])
    expect_true(all(pcx[is.spike,]==0))
    
    X3 <- denoisePCA(X, technical=fit$trend, value="lowrank", get.spikes=TRUE)
    ref <- denoisePCA(exprs(X), technical=fit$trend, value="lowrank")
    pcx <- assay(X3, "lowrank")
    expect_equal(pcx, ref)
})
