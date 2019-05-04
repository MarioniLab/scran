# This checks the denoisePCA function.
# require(scran); require(testthat); source("setup.R"); source("test-parallel-pca.R")

set.seed(1000)
ngenes <- 1000
ncells <- 200
npops <- 10

mu <- matrix(rnorm(npops*ngenes), ncol=npops)
lcounts <- mu[,sample(ncol(mu), ncells, replace=TRUE)]
lcounts <- lcounts + rnorm(length(lcounts))

set.seed(1001)
test_that("parallelPCA works as expected", {
    threshold <- 0.1
    pcs <- parallelPCA(lcounts, value="pca", niters=20, threshold=threshold)
    permuted <- attr(pcs, "permuted.percentVar")
    original <- attr(pcs, "percentVar")

    pvals <- rowMeans(t(permuted) >= original)
    expect_identical(ncol(pcs), min(which(pvals > threshold)-1L))
    expect_true(sum(original) <= 1)
    expect_true(all(rowSums(permuted) <= 1))

    var.exp <- prcomp(t(lcounts))$sdev^2
    frac.exp <- var.exp/sum(var.exp)
    expect_equal(frac.exp[seq_along(original)], original)

    # Respects max.rank when choosing the number of PCs.
    set.seed(100)
    pcs.x <- parallelPCA(lcounts, value="pca", niters=3, max.rank=5, min.rank=0)
    expect_identical(ncol(pcs.x), 5L)
    expect_identical(pcs.x[,], pcs[,1:5])
})

test_that("parallelPCA respects the seed", {
    set.seed(100)
    pcs <- parallelPCA(lcounts, value="pca", niters=3)
    set.seed(100)
    pcs2 <- parallelPCA(lcounts, value="pca", niters=3)
    expect_identical(pcs2, pcs)
    pcs3 <- parallelPCA(lcounts, value="pca", niters=3)
    expect_false(identical(pcs3, pcs))

    # With irlba:
    set.seed(100)
    ipcs <- parallelPCA(lcounts, value="pca", niters=3, BSPARAM=BiocSingular::IrlbaParam(), max.rank=10)
    set.seed(100)
    ipcs2 <- parallelPCA(lcounts, value="pca", niters=3, BSPARAM=BiocSingular::IrlbaParam(), max.rank=10)
    expect_identical(pcs, pcs2)
    expect_identical(ncol(ipcs), ncol(pcs))

    # With parallelization.
    BPPARAM <- safeBPParam(3) # define BEFORE set.seed, otherwise this sets its own seed.
    set.seed(100)
    alt <- parallelPCA(lcounts, value="pca", niters=3, BPPARAM=BPPARAM)
    expect_identical(alt, pcs)
})

set.seed(1002)
test_that("parallelPCA's C++ code works as expected", {
    trans <- t(lcounts)
    shuffled <- .Call(scran:::cxx_shuffle_matrix, trans, 1, 1L)
    expect_false(identical(shuffled, trans))
    expect_identical(apply(shuffled, 2, sort), apply(trans, 2, sort))

    # Is reproducible.
    shuffled2 <- .Call(scran:::cxx_shuffle_matrix, trans, 1, 1L)
    expect_identical(shuffled, shuffled2)

    # Responds to the seed.
    shuffled3 <- .Call(scran:::cxx_shuffle_matrix, trans, 2, 1L)
    expect_false(identical(shuffled, shuffled3))
    expect_identical(apply(shuffled, 2, sort), apply(shuffled3, 2, sort))

    # Responds to the stream.
    shuffled4 <- .Call(scran:::cxx_shuffle_matrix, trans, 1, 2L)
    expect_false(identical(shuffled, shuffled4))
    expect_identical(apply(shuffled, 2, sort), apply(shuffled4, 2, sort))

    # Matches the result from the test shuffler.
    for (i in 1:10) {
        vals <- cbind(runif(i*10))
        seed <- i*100
        stream <- i
        expect_identical(
            .Call(scran:::cxx_shuffle_matrix, vals, seed, stream),
            scramble_matrix(vals, seed=seed, stream=stream)
        )
    }
})

set.seed(1003)
test_that("parallelPCA works with different settings", {
    # Responds correctly to subsetting.
    chosen <- sample(nrow(lcounts), 100)
    set.seed(100)
    pcs <- parallelPCA(lcounts, value="pca", niters=3, subset.row=chosen)
    set.seed(100)
    pcs2 <- parallelPCA(lcounts[chosen,,drop=FALSE], value="pca", niters=3)
    expect_identical(pcs, pcs2)

    # Checking that subsetting and projections work.
    set.seed(100)
    lr <- parallelPCA(lcounts, value="lowrank", niters=3)

    combined <- rbind(lcounts, lcounts[1:10,])
    chosen <- seq_len(nrow(lcounts))
    set.seed(100)
    lr3 <- parallelPCA(combined, value="lowrank", niters=3, subset.row=chosen)
    expect_equal(lr[,], lr3[chosen,])
    expect_equal(lr3[1:10,], lr3[1:10+nrow(lcounts),])
})

test_that("parallelPCA works with SingleCellExperiment inputs", {
    X <- SingleCellExperiment(list(logcounts=lcounts))
    set.seed(120)
    X2 <- parallelPCA(X, niters=3)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL

    set.seed(120)
    pcs <- parallelPCA(lcounts, niters=3)
    expect_identical(pcs, pcx)
    
    # Adding spike-ins, but setting get.spikes=TRUE to avoid subsetting.
    is.spike <- 1:20
    isSpike(X, "Spike") <- is.spike
    set.seed(120)
    X2 <- parallelPCA(X, niters=3, get.spikes=TRUE)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL
    expect_identical(pcx, pcs)
    
    # Trying again, without protecting against spike-ins.
    set.seed(120)
    X2 <- parallelPCA(X, niters=3)
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL

    set.seed(120)
    not.spike <- setdiff(seq_len(ngenes), is.spike)
    pcs <- parallelPCA(lcounts, niters=3, subset.row=not.spike)
    expect_identical(pcs, pcx)

    # Checking lowrank calculations.
    set.seed(120)
    X3 <- parallelPCA(X, niters=3, value="lowrank")

    set.seed(120)
    ref <- parallelPCA(exprs(X)[not.spike,], niters=3, value="lowrank")
    pcx <- assay(X3, "lowrank")
    expect_equal(pcx[not.spike,], ref[,])
    expect_true(all(pcx[is.spike,]==0))
    
    # ... and with all spikes.
    set.seed(120)
    X3 <- parallelPCA(X, niters=3, value="lowrank", get.spikes=TRUE)

    set.seed(120)
    ref <- parallelPCA(exprs(X), niters=3, value="lowrank")
    pcx <- assay(X3, "lowrank")
    expect_equal(pcx, ref)
})
