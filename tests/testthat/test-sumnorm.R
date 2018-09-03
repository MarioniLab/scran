# This tests out the methods related to the deconvolution method.
# require(scran); require(testthat); source("test-sumnorm.R")

ncells <- 200
ngenes <- 1000

set.seed(20000)
test_that("computeSumFactors work correctly on trivial examples", {
    count.sizes <- rnbinom(ncells, mu=100, size=5)
    dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
    out <- computeSumFactors(dummy)
    expect_equal(out, count.sizes/mean(count.sizes))

    # Adding some DE genes.
    count.sizes <- rnbinom(ncells, mu=100, size=5)
    dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
    is.de <- sample(ngenes, 100)
    dummy[is.de,] <- rnbinom(ncells*length(is.de), mu=100, size=1)
    out <- computeSumFactors(dummy)
    expect_equal(out, count.sizes/mean(count.sizes))

    count.sizes <- rnbinom(ncells, mu=100, size=5)
    dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
    is.de <- sample(ngenes, 400)
    dummy[is.de,] <- rnbinom(ncells*length(is.de), mu=100, size=1)
    out <- computeSumFactors(dummy)
    expect_equal(out, count.sizes/mean(count.sizes))
})

set.seed(20001)
test_that("ring construction is correct", {
    lib.sizes <- runif(100)
    out <- scran:::.generateSphere(lib.sizes)
    r <- rank(lib.sizes)
    expect_identical(r[out][1:50], 1:50*2-1L) # All odd ranks
    expect_identical(r[out][51:100], 50:1*2) # All even ranks
    expect_identical(r[out][1:100], r[out][101:200]) # Repeated for easy windowing

    lib.sizes <- runif(101)
    out <- scran:::.generateSphere(lib.sizes)
    r <- rank(lib.sizes)
    expect_identical(r[out][1:51], 1:51*2-1L) # All odd ranks
    expect_identical(r[out][52:101], 50:1*2) # All even ranks
    expect_identical(r[out][1:101], r[out][102:202]) # Repeated for easy windowing
})

coreCheck <- function(x, sphere, pool.sizes)
# Mocking up the core function for creating the linear system.
{
    ave.cell <- scater::calcAverage(x)

    #  Manually running through these.
    all.mat <- all.vals <- vector("list", sum(pool.sizes)*ncells)
    i <- 1L
    for (s in pool.sizes) {
        for (w in seq_len(ncells)) {
            chosen <- sphere[w+seq_len(s)-1L]

            current <- integer(ncells)
            current[chosen] <- 1L
            all.mat[[i]] <- current

            ratios <- rowSums(x[,chosen,drop=FALSE])/ave.cell
            all.vals[[i]] <- median(ratios)
            i <- i+1L
        }
    }

    # Adding the low weight additional equations.
    extra.mat <- diag(ncells)*sqrt(scran:::LOWWEIGHT)
    extra.val <- apply(x/ave.cell, 2, median)*sqrt(scran:::LOWWEIGHT)
    final.mat <- rbind(do.call(rbind, all.mat), extra.mat)
    final.val <- c(unlist(all.vals), extra.val)

    # Checking equality with the core output.
    core <- scran:::.create_linear_system(x, ave.cell=ave.cell, sphere=sphere, pool.sizes=pool.sizes)
    last.set <- seq_len(ncells) + nrow(core$design) - ncells # need some reordering for the last set
    core$design[last.set,][sphere[seq_len(ncells)],] <- core$design[last.set,]
    core$output[last.set][sphere[seq_len(ncells)]] <- core$output[last.set]

    final.mat <- as(final.mat, "dgCMatrix")
    expect_equal(final.mat, core$design)
    expect_equal(final.val, core$output)
}

set.seed(20003)
test_that("construction of the linear system agrees with a reference implementation", {
    pool.sizes <- seq(20, 100, 5)
    ngenes <- 250

    x <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    sphere <- scran:::.generateSphere(runif(ncells))
    coreCheck(x, sphere=sphere, pool.sizes=pool.sizes)

    # Repeating with even numbers of genes and no ties to check the median calculation.
    x <- matrix(rgamma(ngenes*ncells, 10, 10), nrow=ngenes, ncol=ncells)
    sphere <- scran:::.generateSphere(runif(ncells))
    coreCheck(x, sphere=sphere, pool.sizes=pool.sizes)

    # Repeating with odd numbers of genes and no ties.
    x <- matrix(rgamma((ngenes+1)*ncells, 10, 10), nrow=ngenes+1, ncol=ncells)
    sphere <- scran:::.generateSphere(runif(ncells))
    coreCheck(x, sphere=sphere, pool.sizes=pool.sizes)
})

####################################################################################################

library(Matrix)
sumInR <- function(x, sizes, center=TRUE, min.mean=0)
# Creating a quick R implementation for comparison.
{
    lib.sizes <- colSums(x)
    ref <- scater::calcAverage(x)
    keep <- ref >= pmax(1e-8, min.mean)
    x <- x[keep,,drop=FALSE]
    ref <- ref[keep]

    ncells <- length(lib.sizes)
    o <- scran:::.generateSphere(lib.sizes)
    all.mat <- all.vals <- vector("list", sum(sizes)*ncells)

    # Running the core function directly.
    core <- scran:::.create_linear_system(x, ave.cell=ref, sphere=o, pool.sizes=sizes)
    final.mat <- core$design
    final.val <- core$output

    nf <- Matrix::solve(Matrix::qr(final.mat), final.val)
    nf <- as.numeric(nf)
    if (center) {
        nf <- nf/mean(nf)
    }
    return(nf)
}

set.seed(20004)
test_that("computeSumFactors agrees with a reference implementation", {
    ngenes2 <- 200
    x <- matrix(rpois(ngenes2*ncells, lambda=10), nrow=ngenes2, ncol=ncells)
    sizes <- seq(20, 100, 5)
    ref <- sumInR(x, sizes)
    obs <- computeSumFactors(x, sizes=sizes, min.mean=0)
    expect_equal(ref, obs)

    # Works if we throw in some zeroes throughout, to test the default filtering.
    x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
    x[sample(nrow(x), 100),] <- 0L
    ref <- sumInR(x, sizes)
    obs <- computeSumFactors(x, sizes=sizes, min.mean=0)
    expect_equal(ref, obs)

    # Works with subsetting.
    x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
    subset.row <- sample(nrow(x), 100)
    obs <- computeSumFactors(x, subset.row=subset.row, sizes=sizes, min.mean=0)
    ref <- computeSumFactors(x[subset.row,], sizes=sizes, min.mean=0)
    expect_equal(ref, obs)
})

set.seed(20005)
test_that("computeSumFactors correctly ignores low-abundance genes", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=1), nrow=ngenes, ncol=ncells)
    sizes <- seq(20, 100, 5)

    # Can't subset 'dummy' directly for testing, as that would change the library sizes.
    outA <- computeSumFactors(dummy, min.mean=0.5, sizes=sizes)
    expect_equal(outA, sumInR(dummy, sizes=sizes, min.mean=0.5))
    outB <- computeSumFactors(dummy, min.mean=1, sizes=sizes)
    expect_equal(outB, sumInR(dummy, sizes=sizes, min.mean=1))

    expect_false(isTRUE(all.equal(outA, outB))) # ensure it's not trivial equality.

    # Interacts properly with the subsetting.
    out <- computeSumFactors(dummy, min.mean=1, subset.row=1:500, sizes=sizes)
    expect_equal(out, sumInR(dummy[1:500,], sizes=sizes, min.mean=1))
})

####################################################################################################

set.seed(20006)
test_that("computeSumFactors behaves correctly with clustering", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    clusters <- rep(1:2, 100)
    sizes <- seq(20, 100, 5)

    obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters, min.mean=0, ref.clust=1)
    ref1 <- sumInR(dummy[,clusters==1], sizes, center=FALSE) # Avoid centering, as this destroys relative information.
    ref2 <- sumInR(dummy[,clusters==2], sizes, center=FALSE)

    pseudo1 <- scater::calcAverage(dummy[,clusters==1])
    pseudo2 <- scater::calcAverage(dummy[,clusters==2])
    ref2 <- ref2*median(pseudo2/pseudo1)

    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

set.seed(200061)
test_that("computeSumFactors behaves correctly with clustering and a mean threshold", {
    ldummy <- matrix(rpois(ncells*ngenes, lambda=1), nrow=ngenes, ncol=ncells)
    clusters <- rep(1:2, 100)
    sizes <- seq(20, 100, 5)

    l1 <- ldummy[,clusters==1]
    l2 <- ldummy[,clusters==2]

    obs <- computeSumFactors(ldummy, sizes=sizes, cluster=clusters, min.mean=1, ref.clust=1)
    ref1 <- sumInR(l1, sizes, center=FALSE, min.mean=1)
    ref2 <- sumInR(l2, sizes, center=FALSE, min.mean=1)

    ave1 <- scater::calcAverage(l1)
    ave2 <- scater::calcAverage(l2)
    grand <- scater::calcAverage(cbind(ave1, ave2))
    expect_equal(grand, (ave1/sum(ave1) + ave2/sum(ave2))/2 * (sum(ave1) + sum(ave2))/2) # check calculation.

    keep <- grand >= 1 # The grand mean applies during re-scaling across clusters.
    ref2 <- ref2*median(ave2[keep]/ave1[keep])

    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

set.seed(20007)
test_that("computeSumFactors correctly subsets 'sizes' for small clusters", {
    # Trying with not-quite-enough cells in one cluster.
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    clusters <- rep(1:2, c(80, 120))
    sizes <- seq(20, 100, 5)

    expect_warning(obs <- computeSumFactors(dummy[,clusters==1], sizes=sizes, min.mean=0), "not enough cells in at least one cluster")
    ref1 <- sumInR(dummy[,clusters==1], sizes[sizes<=sum(clusters==1)], center=FALSE)
    expect_equal(ref1/mean(ref1), obs)

    expect_warning(obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters, min.mean=0), "not enough cells in at least one cluster")
    ref2 <- sumInR(dummy[,clusters==2], sizes, center=FALSE) # Ensure that second cluster isn't affected by subsetting of sizes.
    pseudo1 <- scater::calcAverage(dummy[,clusters==1])
    pseudo2 <- scater::calcAverage(dummy[,clusters==2])
    ref2 <- ref2 * median(pseudo2/pseudo1)

    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

set.seed(20008)
test_that("computeSumFactors correctly limits cluster sizes", {
    # Checking that it does the job inside the function.
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    obs <- computeSumFactors(dummy, max.cluster.size=100)
    expect_equal(obs, computeSumFactors(dummy, clusters=rep(1:2, length.out=ncol(dummy))))

    # Checking that the size-capping function works.
    clusters <- sample(1:5, 51, p=1:5, replace=TRUE)
    out <- scran:::.limit_cluster_size(clusters, 10)
    expect_true(all(table(out) <= 10L))
    expect_false(identical(out, clusters))
    expect_true(length(unique(paste0(clusters, out)))==length(unique(out))) # nested

    # Checking that it works with factors.
    clusters <- factor(integer(100))
    out <- scran:::.limit_cluster_size(clusters, 6)
    expect_true(all(table(out) <= 6L))
    expect_false(identical(out, clusters))
    expect_true(length(unique(paste0(clusters, out)))==length(unique(out))) # nested
})

set.seed(20009)
test_that("computeSumFactors is correct with clustering in majority-DE cases", {
    ncells <- 600
    ngenes <- 200
    count.sizes <- rnbinom(ncells, mu=100, size=5)
    multiplier <- seq_len(ngenes)/100
    dummy <- outer(multiplier, count.sizes)

    # Most genes (120 out of 200) are DE in at least one cluster.
    known.clusters <- sample(3, ncells, replace=TRUE)
    dummy[1:40,known.clusters==1L] <- 0
    dummy[41:80,known.clusters==2L] <- 0
    dummy[81:120,known.clusters==3L] <- 0

    out <- computeSumFactors(dummy, cluster=known.clusters)
    expect_equal(out, count.sizes/mean(count.sizes)) # Even though there is a majority of DE, each pair of clusters is still okay.

    out1 <- computeSumFactors(dummy, cluster=known.clusters, ref=1)
    expect_equal(out, out1)
    out2 <- computeSumFactors(dummy, cluster=known.clusters, ref=2)
    expect_equal(out, out2)
    out3 <- computeSumFactors(dummy, cluster=known.clusters, ref=3)
    expect_equal(out, out3)

    expect_error(computeSumFactors(dummy, cluster=known.clusters, ref=0), "'ref.clust' value not in 'clusters'")
})

####################################################################################################

set.seed(20010)
test_that("computeSumFactors works on SingleCellExperiment objects", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    rownames(dummy) <- paste0("X", seq_len(ngenes))
    X <- SingleCellExperiment(list(counts=dummy))
    out <- computeSumFactors(X)
    expect_equal(unname(sizeFactors(out)), computeSumFactors(dummy))

    # Trying with and without spike-ins.
    is.spike <- sample(ngenes, 100)
    isSpike(X, "Spike") <- is.spike
    out2 <- computeSumFactors(X)
    expect_equal(unname(sizeFactors(out2)), computeSumFactors(dummy[-is.spike,]))

    out3 <- computeSumFactors(X, get.spikes=TRUE)
    expect_equal(sizeFactors(out3), sizeFactors(out))

    # Combined with subsetting.
    out4 <- computeSumFactors(X, get.spikes=FALSE, subset.row=1:500)
    expect_equal(unname(sizeFactors(out4)), computeSumFactors(dummy[setdiff(1:500, is.spike),]))

    out5 <- computeSumFactors(X, get.spikes=TRUE, subset.row=1:500)
    expect_equal(unname(sizeFactors(out5)), computeSumFactors(dummy[1:500,]))
})

set.seed(20011)
test_that("other solving options work properly", {
    dummy <- matrix(rpois(ncells*ngenes, lambda=10), nrow=ngenes, ncol=ncells)
    out <- computeSumFactors(dummy)
    if (.Platform$OS.type!="windows") { # Because limSolve doesn't build on Windows, apparently.
        outx <- computeSumFactors(dummy, positive=TRUE)
        expect_true(all(abs(outx -  out) < 1e-3)) # need to be a bit generous here, the solution code is different.
    }
    expect_warning(outx <- computeSumFactors(dummy, errors=TRUE), "errors=TRUE is no longer supported")
    expect_equal(as.numeric(outx), out)
})

set.seed(200111)
test_that("computeSumFactors works properly on alternative representations", {
    library(Matrix)
    X <- as(matrix(rpois(100000, lambda=1), ncol=100), "dgCMatrix")
    X_ <- as.matrix(X)
    
    library(HDF5Array)
    Y <- as(matrix(rpois(100000, lambda=5), ncol=100), "HDF5Array")
    Y_ <- as.matrix(Y)

    sf1 <- computeSumFactors(X_, min.mean=0)
    sf2 <- computeSumFactors(X, min.mean=0)
    expect_identical(sf1, sf2)
    
    sf1 <- computeSumFactors(Y_, min.mean=0)
    sf2 <- computeSumFactors(Y, min.mean=0)
    expect_identical(sf1, sf2)
})

set.seed(20012)
test_that("computeSumFactors throws errors correctly", {
    dummy <- matrix(rpois(ncells*ngenes, lambda=10), nrow=ngenes, ncol=ncells)
    expect_error(computeSumFactors(dummy, min.mean=NULL), "turn off abundance filtering")
    expect_error(computeSumFactors(dummy[,0,drop=FALSE]), "not enough cells in at least one cluster")
    expect_error(computeSumFactors(dummy[0,,drop=FALSE]), "insufficient features for median calculations")
    expect_error(computeSumFactors(dummy, sizes=c(10, 10, 20)), "'sizes' are not unique")
    expect_error(computeSumFactors(dummy, clusters=integer(0)), "'ncol(x)' is not equal to 'length(clusters)'", fixed=TRUE)
})
