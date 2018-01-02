# This tests out the normalization methods in scran - specifically, compute*Factors and normalize().
# require(scran); require(testthat); source("test-normalize.R")

set.seed(20000)
ncells <- 200
ngenes <- 1000

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

# Checking the ring construction.

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

# Checking the subset and divide function is working properly.

test_that("subset and division is correct", {
    x <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    subset.row <- sample(ngenes, 500)
    subset.col <- sample(ncells, 100)
    cur.out <- .Call(scran:::cxx_subset_and_divide, x, subset.row-1L, subset.col-1L) 

    chosen <- x[subset.row,subset.col]
    expect_equal(cur.out[[1]], colSums(chosen))
    divved <- t(t(chosen)/colSums(chosen))
    expect_equal(cur.out[[2]], divved)
    expect_equal(cur.out[[3]], rowMeans(divved))

    # Checking the logic of the mean abundance.
    expect_equal(cur.out[[3]]*mean(cur.out[[1]]), scater::calcAverage(chosen))
})

# Checking the core function for creating the linear system.

coreCheck <- function(x, sphere, pool.sizes) {
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

test_that("construction of the linear system agrees with a reference implementation", {
    pool.sizes <- seq(20, 100, 5)
    ngenes <- 250

    set.seed(3000)
    x <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    sphere <- scran:::.generateSphere(runif(ncells))
    coreCheck(x, sphere=sphere, pool.sizes=pool.sizes)

    # Repeating with even numbers of genes and no ties to check the median calculation.
    set.seed(3001)
    x <- matrix(rgamma(ngenes*ncells, 10, 10), nrow=ngenes, ncol=ncells)
    sphere <- scran:::.generateSphere(runif(ncells))
    coreCheck(x, sphere=sphere, pool.sizes=pool.sizes)
    
    # Repeating with odd numbers of genes and no ties.
    set.seed(3002)
    x <- matrix(rgamma((ngenes+1)*ncells, 10, 10), nrow=ngenes+1, ncol=ncells)
    sphere <- scran:::.generateSphere(runif(ncells))
    coreCheck(x, sphere=sphere, pool.sizes=pool.sizes)
})

####################################################################################################

# Creating a quick R implementation for comparison.

library(Matrix)
sumInR <- function(x, sizes, center=TRUE, min.mean=0) {
    keep <- scater::calcAverage(x) >= pmax(1e-8, min.mean)
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
    ref <- sumInR(x[subset.row,,drop=FALSE], sizes)
    obs <- computeSumFactors(x, subset.row=subset.row, sizes=sizes, min.mean=0)
    expect_equal(ref, obs)
})

####################################################################################################

# Trying it out with other options.

dummy <- matrix(rpois(ncells*ngenes, lambda=10), nrow=ngenes, ncol=ncells)
test_that("other solving options work properly", {
    out <- computeSumFactors(dummy)
    if (.Platform$OS.type!="windows") { # Because limSolve doesn't build on Windows, apparently.
        outx <- computeSumFactors(dummy, positive=TRUE)
        expect_true(all(abs(outx -  out) < 1e-3)) # need to be a bit generous here, the solution code is different.
    }
    expect_warning(outx <- computeSumFactors(dummy, errors=TRUE), "errors=TRUE is no longer supported")
    expect_equal(as.numeric(outx), out)
})

# Checking the the clustering works as expected.

test_that("computeSumFactors behaves correctly with clustering", {
    clusters <- rep(1:2, 100)
    sizes <- seq(20, 100, 5)
    obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters, min.mean=0, ref.clust=1)
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

    # Checking that we get the same performance with a mean threshold.
    ldummy <- matrix(rpois(ncells*ngenes, lambda=1), nrow=ngenes, ncol=ncells)
    l1 <- ldummy[,clusters==1]
    l2 <- ldummy[,clusters==2]

    obs <- computeSumFactors(ldummy, sizes=sizes, cluster=clusters, min.mean=1, ref.clust=1)
    ref1 <- sumInR(l1, sizes, center=FALSE, min.mean=1) 
    ref2 <- sumInR(l2, sizes, center=FALSE, min.mean=1)

    adj1 <- t(t(l1)/colSums(l1))
    adj2 <- t(t(l2)/colSums(l2))
    pseudo1 <- rowMeans(adj1)
    pseudo2 <- rowMeans(adj2)

    ave1 <- scater::calcAverage(l1)
    ave2 <- scater::calcAverage(l2)
    keep <- (ave1+ave2)/2 >= 1 # The grand mean applies during re-scaling across clusters.
    ref2 <- ref2*median(pseudo2[keep]/pseudo1[keep])
    
    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

test_that("computeSumFactors correctly subsets 'sizes' for small clusters", {
    # Trying with not-quite-enough cells in one cluster.
    clusters <- rep(1:2, c(80, 120))
    sizes <- seq(20, 100, 5)
    expect_warning(obs <- computeSumFactors(dummy[,clusters==1], sizes=sizes, min.mean=0), "not enough cells in at least one cluster")
    ref1 <- sumInR(dummy[,clusters==1], sizes[sizes<=sum(clusters==1)], center=FALSE) 
    expect_equal(ref1/mean(ref1), obs)
    
    expect_warning(obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters, min.mean=0), "not enough cells in at least one cluster")
    ref2 <- sumInR(dummy[,clusters==2], sizes, center=FALSE) # Ensure that second cluster isn't affected by subsetting of sizes.
    adj <- t(t(dummy)/colSums(dummy))
    pseudo1 <- rowMeans(adj[,clusters==1])
    pseudo2 <- rowMeans(adj[,clusters==2])
    ref2 <- ref2 * median(pseudo2/pseudo1)
    
    ref <- numeric(ncells)
    ref[clusters==1] <- ref1
    ref[clusters==2] <- ref2
    ref <- ref/mean(ref)
    expect_equal(ref, obs)
})

set.seed(20004)
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

# Trying it out on a SingleCellExperiment object.

set.seed(20001)
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

# Testing that the minimum mean specifications work.

set.seed(20002)
test_that("computeSumFactors correctly detects low-abundance genes", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=1), nrow=ngenes, ncol=ncells)
    sizes <- seq(20, 100, 5)
    
    # Can't subset 'dummy' directly for testing, as that would change the library sizes.
    out <- computeSumFactors(dummy, min.mean=0.5)
    expect_equal(out, sumInR(dummy, sizes=sizes, min.mean=0.5))
    out <- computeSumFactors(dummy, min.mean=1)
    expect_equal(out, sumInR(dummy, sizes=sizes, min.mean=1))

    # Interacts properly with the subsetting.
    out <- computeSumFactors(dummy, min.mean=1, subset.row=1:500)
    expect_equal(out, sumInR(dummy[1:500,], sizes=sizes, min.mean=1))
})

# Throwing in some silly inputs.

test_that("computeSumFactors throws errors correctly", {
    expect_error(computeSumFactors(dummy[,0,drop=FALSE]), "not enough cells in at least one cluster")
    expect_error(computeSumFactors(dummy[0,,drop=FALSE]), "cells should have non-zero library sizes")
    expect_error(computeSumFactors(dummy, sizes=c(10, 10, 20)), "'sizes' are not unique")
    expect_error(computeSumFactors(dummy, clusters=integer(0)), "'x' ncols is not equal to 'clusters' length")
})

####################################################################################################

# Checking out the behaviour of the computeSpikeFactors function.

set.seed(20003)
ncells <- 200
ngenes <- 1000

dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
is.spike <- rbinom(ngenes, 1, 0.7)==0L
dummy[is.spike,] <- matrix(rnbinom(sum(is.spike)*ncells, mu=20, size=5), ncol=ncells, nrow=sum(is.spike), byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
isSpike(X, "MySpike") <- is.spike

test_that("computeSpikeFactors calculates spike-based size factors correctly", {
    out <- computeSpikeFactors(X)
    ref <- colSums(dummy[is.spike,])
    expect_equal(unname(sizeFactors(out)), ref/mean(ref))
    expect_equal(sizeFactors(out), sizeFactors(out, type="MySpike"))
})

test_that("computeSpikeFactors works with multiple spike-in sets", {
    subset <- split(which(is.spike), rep(1:2, length.out=sum(is.spike)))
    isSpike(X, "MySpike") <- subset[[1]]
    isSpike(X, "SecondSpike") <- subset[[2]]
    out <- computeSpikeFactors(X)

    out.sub <- computeSpikeFactors(X, type="MySpike") # Sanity check, to make sure that it's calculating it differently for each spike-in.
    subref <- colSums(dummy[subset[[1]],])
    expect_equal(unname(sizeFactors(out.sub)), subref/mean(subref))
    expect_equal(sizeFactors(out.sub), sizeFactors(out.sub, type="MySpike"))
    expect_identical(sizeFactors(out.sub, type="SecondSpike"), NULL)
    
    out2 <- computeSpikeFactors(X)
    expect_equal(sizeFactors(out), sizeFactors(out2))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))
    
    out2 <- computeSpikeFactors(X, type=c("MySpike", "SecondSpike"))
    expect_equal(sizeFactors(out), sizeFactors(out2))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
    expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))
})

test_that("computeSpikeFactors responds correctly to general.use", {
    sizeFactors(X) <- 1
    out <- computeSpikeFactors(X, general.use=FALSE)
    expect_equal(unname(sizeFactors(out)), rep(1, ncells))
    ref <- colSums(dummy[is.spike,])
    expect_equal(unname(sizeFactors(out, type="MySpike")), ref/mean(ref))
})

test_that("computeSpikeFactors fails correctly on silly inputs", {
    expect_error(out <- computeSpikeFactors(X[0,]), "no spike-in transcripts present in 'x'")

    alt.X <- X
    counts(alt.X)[] <- 0L
    expect_warning(out <- computeSpikeFactors(alt.X), "zero spike-in counts during spike-in normalization")
    expect_identical(unname(sizeFactors(out)), rep(NaN, ncol(out)))

    # Checking that it correctly returns nothing.
    out <- computeSpikeFactors(X[,0])
    expect_identical(unname(sizeFactors(out)), numeric(0))
})

####################################################################################################

# Checking out the behaviour of the normalize() function.

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))

X <- SingleCellExperiment(list(counts=dummy))
ref <- colSums(dummy)
sizeFactors(X) <- ref

test_that("normalize() from scater works on endogenous genes", {
    out <- normalize(X)
    sf <- ref/mean(ref)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1))
    out <- normalize(X, log_exprs_offset=3)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+3))
    
    # Repeating with different set of size factors.
    ref <- runif(ncells, 10, 20)
    sizeFactors(X) <- ref
    out <- normalize(X)
    sf <- ref/mean(ref)
    expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1)) 

    # Checking that size factor centering works properly.    
    expect_equivalent(sf, sizeFactors(out))
    Xb <- X
    sizeFactors(Xb) <- ref
    outb <- normalize(Xb, centre_size_factors=FALSE)
    expect_equivalent(ref, sizeFactors(outb))
    expect_equivalent(exprs(out), exprs(outb))
})

# Now adding some controls.

test_that("normalize from scater works on spike-in genes", {
    out <- normalize(X)
    chosen <- rbinom(ngenes, 1, 0.7)==0L
    isSpike(X, "whee") <- chosen
    expect_warning(X3 <- normalize(X), "spike-in set 'whee'")
    expect_equal(exprs(out), exprs(X3))

    # Checking that it correctly uses the spike-in size factors.    
    sizeFactors(X, type="whee") <- colSums(counts(X)[chosen,])
    expect_warning(X4 <- normalize(X), NA) # i.e., no warning.
    expect_equivalent(exprs(out)[!chosen,], exprs(X4)[!chosen,])
    ref <- sizeFactors(X, type="whee")
    sf <- ref/mean(ref)
    expect_equivalent(exprs(X4)[chosen,], log2(t(t(dummy[chosen,])/sf)+1))
    
    expect_equivalent(sizeFactors(X4, type="whee"), sf)
    X4b <- normalize(X, centre_size_factors=FALSE)
    expect_equivalent(sizeFactors(X4b, type="whee"), sizeFactors(X, type="whee"))
    expect_equivalent(exprs(X4), exprs(X4b))
})

# Checking out silly inputs.

expect_equal(unname(dim(normalize(X[,0,drop=FALSE]))), c(ngenes, 0L))
expect_equal(unname(dim(normalize(X[0,,drop=FALSE]))), c(0L, ncells)) 


