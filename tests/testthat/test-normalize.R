# This tests out the normalization methods in scran - specifically, compute*Factors and normalize().
# require(scran); require(testthat); source("test-normalize.R")

set.seed(20000)
ncells <- 200
ngenes <- 1000
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

# Checking the ring construction.

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

# Creating an R-only implementation for comparison.

sumInR <- function(x, sizes, center=TRUE) {
    lib.sizes <- colSums(x)
    x <- t(t(x)/lib.sizes)
    ref <- rowMeans(x)

    keep <- ref > 0
    ref <- ref[keep]
    x <- x[keep,,drop=FALSE]

    ncells <- length(lib.sizes)
    o <- scran:::.generateSphere(lib.sizes)
    all.mat <- all.vals <- vector("list", sum(sizes)*ncells) 
    i <- 1L

    for (s in sizes) {
        for (w in seq_len(ncells)) {
            chosen <- o[w+seq_len(s)-1L]

            current <- integer(ncells)
            current[chosen] <- 1L
            all.mat[[i]] <- current

            ratios <- rowSums(x[,chosen,drop=FALSE])/ref
            all.vals[[i]] <- median(ratios)
            i <- i+1L
        }
    }

    # Adding the low weight additional equations.
    extra.mat <- diag(ncells)*sqrt(scran:::LOWWEIGHT)
    extra.val <- apply(x/ref, 2, median)*sqrt(scran:::LOWWEIGHT)
    final.mat <- rbind(do.call(rbind, all.mat), extra.mat)
    final.val <- c(unlist(all.vals), extra.val)

    nf <- solve(qr(final.mat), final.val)
    sf <- nf * lib.sizes
    if (center) {
        sf <- sf/mean(sf)
    }
    return(sf)
}

ngenes2 <- 200
x <- matrix(rpois(ngenes2*ncells, lambda=10), nrow=ngenes2, ncol=ncells)
sizes <- seq(20, 100, 5)
ref <- sumInR(x, sizes)
obs <- computeSumFactors(x, sizes=sizes)
expect_equal(ref, obs)

x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
x[sample(nrow(x), 100),] <- 0L # Throwing in some zeroes.
ref <- sumInR(x, sizes)
obs <- computeSumFactors(x, sizes=sizes, mean.warn=FALSE) # shutting up the warnings here.
expect_equal(ref, obs)

x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
subset.row <- sample(nrow(x), 100)
ref <- sumInR(x[subset.row,,drop=FALSE], sizes)
obs <- computeSumFactors(x, subset.row=subset.row, sizes=sizes)
expect_equal(ref, obs)

####################################################################################################

# Trying it out with other options.

dummy <- matrix(rpois(ncells*ngenes, lambda=10), nrow=ngenes, ncol=ncells)
out <- computeSumFactors(dummy)
if (.Platform$OS.type!="windows") { # Because limSolve doesn't build on Windows, apparently.
outx <- computeSumFactors(dummy, positive=TRUE)
expect_true(all(abs(outx -  out) < 1e-3)) # need to be a bit generous here, the solution code is different.
}
expect_warning(outx <- computeSumFactors(dummy, errors=TRUE), "errors=TRUE is no longer supported")
expect_equal(as.numeric(outx), out)

# Checking the the clustering works as expected.

clusters <- rep(1:2, 100)
sizes <- seq(20, 100, 5)
obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters)
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

# Trying with not-quite-enough cells in one cluster.

test_that("computeSumFactors correctly subsets 'sizes' for small clusters", {
    clusters <- rep(1:2, c(80, 120))
    sizes <- seq(20, 100, 5)
    expect_warning(obs <- computeSumFactors(dummy[,clusters==1], sizes=sizes), "not enough cells in at least one cluster")
    ref1 <- sumInR(dummy[,clusters==1], sizes[sizes<=80], center=FALSE) 
    expect_equal(ref1/mean(ref1), obs)
    
    expect_warning(obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters), "not enough cells in at least one cluster")
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

# Trying it out on a SCESet object.

set.seed(20001)
test_that("computeSumFactors works on SCESest objects", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=10), nrow=ngenes, ncol=ncells)
    rownames(dummy) <- paste0("X", seq_len(ngenes))
    X <- newSCESet(countData=data.frame(dummy))
    out <- computeSumFactors(X)
    expect_equal(unname(sizeFactors(out)), computeSumFactors(dummy))

    # Trying with and without spike-ins.
    is.spike <- sample(ngenes, 100)
    X <- calculateQCMetrics(X, feature_controls=list(Spike=is.spike))
    setSpike(X) <- "Spike"
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

# Testing that the minimum mean warnings work.

set.seed(20002)
test_that("computeSumFactors correctly detects low-abundance genes", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=0.1), nrow=ngenes, ncol=ncells)
    expect_warning(computeSumFactors(dummy), "low-abundance genes")
    expect_warning(computeSumFactors(dummy, mean.warn=FALSE), NA)

    sf <- colSums(dummy)
    sf <- sf/mean(sf)
    corrected.counts <- t(t(dummy)/sf)
    keep <- rowMeans(corrected.counts) > 0.11 # slightly higher, as library sizes change after filtering.
    expect_warning(computeSumFactors(dummy, subset.row=keep), NA)
})

# Throwing in some silly inputs.

expect_error(computeSumFactors(dummy[,0,drop=FALSE]), "not enough cells in at least one cluster")
expect_error(computeSumFactors(dummy[0,,drop=FALSE]), "cells should have non-zero library sizes")
expect_error(computeSumFactors(dummy, sizes=c(10, 10, 20)), "'sizes' are not unique")
expect_error(computeSumFactors(dummy, clusters=integer(0)), "'x' ncols is not equal to 'clusters' length")

####################################################################################################

# Checking out the behaviour of the computeSpikeFactors function.

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
is.spike <- rbinom(ngenes, 1, 0.7)==0L
dummy[is.spike,] <- matrix(rnbinom(sum(is.spike)*ncells, mu=20, size=5), ncol=ncells, nrow=sum(is.spike), byrow=TRUE)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
X <- calculateQCMetrics(X, list(MySpike=is.spike))
setSpike(X) <- "MySpike"
out <- computeSpikeFactors(X)
ref <- colSums(dummy[is.spike,])
expect_equal(unname(sizeFactors(out)), ref/mean(ref))
expect_equal(sizeFactors(out), sizeFactors(out, type="MySpike"))

# Checking out what happens when you have multiple spike-ins supplied.
X2 <- newSCESet(countData=data.frame(dummy))
subset <- split(which(is.spike), rep(1:2, length.out=sum(is.spike)))
X2 <- calculateQCMetrics(X2, list(MySpike=subset[[1]], SecondSpike=subset[[2]]))
setSpike(X2) <- c("MySpike", "SecondSpike")

out.sub <- computeSpikeFactors(X2, type="MySpike") # Sanity check, to make sure that it's calculating it differently for each spike-in.
subref <- colSums(dummy[subset[[1]],])
expect_equal(unname(sizeFactors(out.sub)), subref/mean(subref))
expect_equal(sizeFactors(out.sub), sizeFactors(out.sub, type="MySpike"))
expect_warning(sizeFactors(out.sub, type="SecondSpike"), "'sizeFactors' have not been set for 'SecondSpike'")

out2 <- computeSpikeFactors(X2)
expect_equal(sizeFactors(out), sizeFactors(out2))
expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))

out2 <- computeSpikeFactors(X2, type=c("MySpike", "SecondSpike"))
expect_equal(sizeFactors(out), sizeFactors(out2))
expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))

# Checking out the general use function.
sizeFactors(X) <- 1
out <- computeSpikeFactors(X, general.use=FALSE)
expect_equal(unname(sizeFactors(out)), rep(1, ncells))
expect_equal(unname(sizeFactors(out, type="MySpike")), ref/mean(ref))

# Breaks if you try to feed it silly inputs.
expect_warning(out <- computeSpikeFactors(X[0,]), "zero spike-in counts during spike-in normalization")
expect_identical(unname(sizeFactors(out)), rep(NaN, ncol(out)))
out <- computeSpikeFactors(X[,0])
expect_identical(unname(sizeFactors(out)), numeric(0))

####################################################################################################

# Checking out the behaviour of the normalize() function.

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))
X <- newSCESet(countData=dummy)

ref <- colSums(dummy)
sizeFactors(X) <- ref
out <- normalize(X)
sf <- ref/mean(ref)
expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1))
out <- normalize(X, logExprsOffset=3)
expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+3))

ref <- runif(ncells, 10, 20)
sizeFactors(X) <- ref
out <- normalize(X)
sf <- ref/mean(ref)
expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1)) 

expect_equivalent(sf, sizeFactors(out))
Xb <- X
sizeFactors(Xb) <- ref
outb <- normalize(Xb, centre_size_factors=FALSE)
expect_equivalent(ref, sizeFactors(outb))
expect_equivalent(exprs(out), exprs(outb))

# Now adding some controls.

chosen <- rbinom(ngenes, 1, 0.7)==0L
X <- calculateQCMetrics(X, feature_controls=list(whee=chosen))
X3 <- normalize(X)
expect_equal(exprs(out), exprs(X3))

Xb <- X
setSpike(Xb) <- "whee"
expect_warning(X3b <- normalize(Xb), "spike-in transcripts in 'whee'")
expect_equal(exprs(X3b), exprs(X3))

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

# Checking out silly inputs.

expect_equal(unname(dim(normalize(X[,0,drop=FALSE]))), c(ngenes, 0L))
expect_equal(unname(dim(normalize(X[0,,drop=FALSE]))), c(0L, ncells)) 


