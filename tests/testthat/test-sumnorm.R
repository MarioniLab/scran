# This tests out the methods related to the deconvolution method.
# require(scran); require(testthat); source("test-sumnorm.R")

ncells <- 200
ngenes <- 1000

set.seed(20000)
test_that("calculateSumFactors work correctly on trivial examples", {
    count.sizes <- rnbinom(ncells, mu=100, size=5)
    dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
    out <- calculateSumFactors(dummy)
    expect_equal(out, count.sizes/mean(count.sizes))

    # Adding some DE genes.
    count.sizes <- rnbinom(ncells, mu=100, size=5)
    dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
    is.de <- sample(ngenes, 100)
    dummy[is.de,] <- rnbinom(ncells*length(is.de), mu=100, size=1)
    out <- calculateSumFactors(dummy)
    expect_equal(out, count.sizes/mean(count.sizes))

    count.sizes <- rnbinom(ncells, mu=100, size=5)
    dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
    is.de <- sample(ngenes, 400)
    dummy[is.de,] <- rnbinom(ncells*length(is.de), mu=100, size=1)
    out <- calculateSumFactors(dummy)
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
    x <- t(t(x)/colSums(x))
    ave.cell <- rowMeans(x)

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
    extra.val <- rep(sqrt(scran:::LOWWEIGHT) / sum(ave.cell), ncells)
    final.mat <- rbind(do.call(rbind, all.mat), extra.mat)
    final.val <- c(unlist(all.vals), extra.val)

    core <- scran:::.create_linear_system(x, ave.cell=ave.cell, sphere=sphere, pool.sizes=pool.sizes)
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
    keep <- scuttle::calculateAverage(x) >= pmax(1e-8, min.mean)
    lib.sizes <- colSums(x)
    x <- t(t(x)/lib.sizes)
    x <- x[keep,,drop=FALSE]
    ref <- rowMeans(x)

    ncells <- length(lib.sizes)
    o <- scran:::.generateSphere(lib.sizes)
    all.mat <- all.vals <- vector("list", sum(sizes)*ncells)

    # Running the core function directly.
    core <- scran:::.create_linear_system(x, ave.cell=ref, sphere=o, pool.sizes=sizes)
    final.mat <- core$design
    final.val <- core$output

    nf <- Matrix::solve(Matrix::qr(final.mat), final.val)
    nf <- as.numeric(nf)
    sf <- nf * lib.sizes
    if (center) {
        sf <- sf/mean(sf)
    }
    return(sf)
}

set.seed(20004)
test_that("calculateSumFactors agrees with a reference implementation", {
    ngenes2 <- 200
    x <- matrix(rpois(ngenes2*ncells, lambda=10), nrow=ngenes2, ncol=ncells)
    sizes <- seq(20, 100, 5)
    ref <- sumInR(x, sizes)
    obs <- calculateSumFactors(x, sizes=sizes, min.mean=0)
    expect_equal(ref, obs)

    # Works if we throw in some zeroes throughout, to test the default filtering.
    x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
    x[sample(nrow(x), 100),] <- 0L
    ref <- sumInR(x, sizes)
    obs <- calculateSumFactors(x, sizes=sizes, min.mean=0)
    expect_equal(ref, obs)

    # Works with subsetting.
    x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
    subset.row <- sample(nrow(x), 100)
    obs <- calculateSumFactors(x, subset.row=subset.row, sizes=sizes, min.mean=0)
    ref <- calculateSumFactors(x[subset.row,], sizes=sizes, min.mean=0)
    expect_equal(ref, obs)
})

set.seed(20005)
test_that("calculateSumFactors correctly ignores low-abundance genes", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=seq_len(ngenes)/ngenes*2), nrow=ngenes, ncol=ncells)
    sizes <- seq(20, 100, 5)

    # Can't subset 'dummy' directly for testing, as that would change the library sizes.
    outA <- calculateSumFactors(dummy, min.mean=0.5, sizes=sizes)
    expect_equal(outA, sumInR(dummy, sizes=sizes, min.mean=0.5))
    outB <- calculateSumFactors(dummy, min.mean=1, sizes=sizes)
    expect_equal(outB, sumInR(dummy, sizes=sizes, min.mean=1))

    expect_false(isTRUE(all.equal(outA, outB))) # ensure it's not trivial equality.
    expect_equal(scuttle::calculateAverage(dummy), colMeans(t(dummy)/colSums(dummy)) * mean(colSums(dummy))) # checking the calculation.

    # Interacts properly with the subsetting.
    out <- calculateSumFactors(dummy, min.mean=1, subset.row=1:500, sizes=sizes)
    expect_equal(out, sumInR(dummy[1:500,], sizes=sizes, min.mean=1))

    # Behaves properly with auto-selection of min.mean.
    expect_identical(calculateSumFactors(dummy, sizes=sizes), calculateSumFactors(dummy, min.mean=0.1, sizes=sizes)) # UMI threshold
    expect_equal(calculateSumFactors(dummy*100, sizes=sizes), calculateSumFactors(dummy, min.mean=1/100, sizes=sizes)) # read threshold
})

set.seed(200051)
test_that("calculateSumFactors responds to scaling requests", {
    truth <- runif(ncells, 0.5, 1.5)
    dummy <- matrix(rpois(ngenes*ncells, lambda=truth), nrow=ngenes, ncol=ncells, byrow=TRUE)

    outA <- calculateSumFactors(dummy, min.mean=0, scaling=NULL)
    outB <- calculateSumFactors(dummy, min.mean=0, scaling=scuttle::librarySizeFactors(dummy))
    expect_equal(outA, outB)

    outC <- calculateSumFactors(dummy, min.mean=0, scaling=truth)
    expect_false(isTRUE(all.equal(outA, outC)))

    # Matching it to what is expected.
    truth.order <- 1 + rank(truth)/1e10 # ensuring the sphere order is the same.
    outD <- calculateSumFactors(t(t(dummy)/(truth/mean(truth))), min.mean=0, scaling=truth.order)
    outD <- outD * truth
    expect_equal(outD/mean(outD), outC/mean(outC))

    # Subsetted properly with clusters.
    clusters <- gl(2, ncells/2)
    outA <- calculateSumFactors(dummy, clusters=clusters, min.mean=0, scaling=NULL)
    outB <- calculateSumFactors(dummy, clusters=clusters, min.mean=0, scaling=scuttle::librarySizeFactors(dummy))
    expect_equal(outA, outB)

    # Throws upon silly inputs.
    expect_error(calculateSumFactors(dummy, min.mean=0, scaling=1), "should be equal")
})

####################################################################################################

set.seed(20006)
test_that("calculateSumFactors behaves correctly with clustering", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    clusters <- rep(1:2, 100)
    sizes <- seq(20, 100, 5)

    obs <- calculateSumFactors(dummy, sizes=sizes, cluster=clusters, min.mean=0, ref.clust=1)
    ref1 <- sumInR(dummy[,clusters==1], sizes, center=FALSE) # Avoid centering, as this destroys relative information.
    ref2 <- sumInR(dummy[,clusters==2], sizes, center=FALSE)

    adj <- t(t(dummy)/colSums(dummy))
    pseudo1 <- rowMeans(adj[,clusters==1])
    pseudo2 <- rowMeans(adj[,clusters==2])
    ref2 <- ref2*median(pseudo2/pseudo1)

    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

set.seed(200061)
test_that("calculateSumFactors behaves correctly with clustering and a mean threshold", {
    ldummy <- matrix(rpois(ncells*ngenes, lambda=1), nrow=ngenes, ncol=ncells)
    clusters <- rep(1:2, 100)
    sizes <- seq(20, 100, 5)

    l1 <- ldummy[,clusters==1]
    l2 <- ldummy[,clusters==2]

    obs <- calculateSumFactors(ldummy, sizes=sizes, cluster=clusters, min.mean=1, ref.clust=1)
    ref1 <- sumInR(l1, sizes, center=FALSE, min.mean=1)
    ref2 <- sumInR(l2, sizes, center=FALSE, min.mean=1)

    adj1 <- t(t(l1)/colSums(l1))
    adj2 <- t(t(l2)/colSums(l2))
    pseudo1 <- rowMeans(adj1)
    pseudo2 <- rowMeans(adj2)

    ave1 <- scuttle::calculateAverage(l1)
    ave2 <- scuttle::calculateAverage(l2)
    grand <- scuttle::calculateAverage(cbind(ave1, ave2))
    expect_equal(grand, (ave1/sum(ave1) + ave2/sum(ave2))/2 * (sum(ave1) + sum(ave2))/2) # check calculation.

    keep <- grand >= 1 # The grand mean applies during re-scaling across clusters.
    ref2 <- ref2*median(pseudo2[keep]/pseudo1[keep])

    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

set.seed(20007)
test_that("calculateSumFactors correctly subsets 'sizes' for small clusters", {
    # Trying with not-quite-enough cells in one cluster.
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    clusters <- rep(1:2, c(80, 120))
    sizes <- seq(20, 100, 5)

    obs <- calculateSumFactors(dummy[,clusters==1], sizes=sizes, min.mean=0)
    ref1 <- sumInR(dummy[,clusters==1], sizes[sizes<=sum(clusters==1)], center=FALSE)
    expect_equal(ref1/mean(ref1), obs)

    # Ensure that second cluster isn't affected by subsetting of sizes.
    obs <- calculateSumFactors(dummy, sizes=sizes, cluster=clusters, min.mean=0)
    ref2 <- sumInR(dummy[,clusters==2], sizes, center=FALSE) 
    adj <- t(t(dummy)/colSums(dummy))
    pseudo1 <- rowMeans(adj[,clusters==1])
    pseudo2 <- rowMeans(adj[,clusters==2])
    ref2 <- ref2 * median(pseudo2/pseudo1)

    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)

    # Degrades to library size normalization.
    subdummy <- dummy[,1:20]
    expect_equal(
            calculateSumFactors(subdummy, sizes=100L, min.mean=0),
            scuttle::librarySizeFactors(subdummy)
    )
})

set.seed(20008)
test_that("calculateSumFactors correctly limits cluster sizes", {
    # Checking that it does the job inside the function.
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    obs <- calculateSumFactors(dummy, max.cluster.size=100)
    expect_equal(obs, calculateSumFactors(dummy, clusters=rep(1:2, length.out=ncol(dummy))))

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
test_that("calculateSumFactors is correct with clustering in majority-DE cases", {
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

    out <- calculateSumFactors(dummy, cluster=known.clusters)
    expect_equal(out, count.sizes/mean(count.sizes)) # Even though there is a majority of DE, each pair of clusters is still okay.

    out1 <- calculateSumFactors(dummy, cluster=known.clusters, ref=1)
    expect_equal(out, out1)
    out2 <- calculateSumFactors(dummy, cluster=known.clusters, ref=2)
    expect_equal(out, out2)
    out3 <- calculateSumFactors(dummy, cluster=known.clusters, ref=3)
    expect_equal(out, out3)

    expect_error(calculateSumFactors(dummy, cluster=known.clusters, ref="0"), "'ref.clust' not in 'clusters'")
})

####################################################################################################

set.seed(20010)
test_that("computeSumFactors works on SingleCellExperiment objects", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    rownames(dummy) <- paste0("X", seq_len(ngenes))
    X <- SingleCellExperiment(list(counts=dummy))
    out <- computeSumFactors(X)
    expect_equal(unname(sizeFactors(out)), calculateSumFactors(dummy))

    # Combined with subsetting.
    out4 <- computeSumFactors(X, subset.row=1:500)
    expect_equal(unname(sizeFactors(out4)), calculateSumFactors(dummy[1:500,]))
})

set.seed(20011)
test_that("setting positive=TRUE behaves properly", {
    lambda <- c(rep(1e-2, 100), 2^rnorm(200))
    dummy <- matrix(rpois(length(lambda)*ngenes, lambda=lambda), nrow=ngenes, ncol=length(lambda), byrow=TRUE)
    expect_warning(out <- calculateSumFactors(dummy), "negative")
    expect_true(all(out > 0)) 

    expect_warning(out2 <- calculateSumFactors(dummy, positive=FALSE), "negative")
    expect_true(any(out2 < 0))

    okay <- out2 > 0
    expect_equal(out[okay]/mean(out[okay]), out2[okay]/mean(out2[okay]))    
})

set.seed(200111)
test_that("calculateSumFactors works properly on alternative representations", {
    library(Matrix)
    X <- as(matrix(rpois(100000, lambda=1), ncol=100), "dgCMatrix")
    X_ <- as.matrix(X)
    
    library(HDF5Array)
    Y <- as(matrix(rpois(100000, lambda=5), ncol=100), "HDF5Array")
    Y_ <- as.matrix(Y)

    sf1 <- calculateSumFactors(X_, min.mean=0)
    sf2 <- calculateSumFactors(X, min.mean=0)
    expect_equal(sf1, sf2)
    
    sf1 <- calculateSumFactors(Y_, min.mean=0)
    sf2 <- calculateSumFactors(Y, min.mean=0)
    expect_equal(sf1, sf2)
})

set.seed(20012)
test_that("calculateSumFactors throws errors correctly", {
    dummy <- matrix(rpois(ncells*ngenes, lambda=10), nrow=ngenes, ncol=ncells)
    expect_error(calculateSumFactors(dummy[,0,drop=FALSE]), "zero cells in one of the clusters")
    expect_error(calculateSumFactors(dummy[0,,drop=FALSE]), "cells should have non-zero library sizes")
    expect_error(calculateSumFactors(dummy, sizes=c(10, 10, 20)), "'sizes' are not unique")
    expect_error(calculateSumFactors(dummy, clusters=integer(0)), "'ncol(x)' is not equal to 'length(clusters)'", fixed=TRUE)
})
