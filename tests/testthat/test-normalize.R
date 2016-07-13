# This tests out the normalization methods in scran - specifically, compute*Factors and normalize().

require(scran); require(testthat)

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

# Trying it out with other options.

outx <- computeSumFactors(dummy, positive=TRUE)
expect_true(all(abs(outx -  out) < 1e-3)) # need to be a bit generous here, the solution code is different.
outx <- computeSumFactors(dummy, errors=TRUE)
expect_equal(as.numeric(outx), out)

expect_identical(names(attributes(outx)), "standard.error") # Checking that the standard errors are equal.
sphere <- scran:::.generateSphere(colSums(dummy))
sizes <- c(20, 40, 60, 80, 100)
use.ave.cell <- rowMeans(dummy)
collected <- scran:::.create_linear_system(dummy, sphere, as.integer(sizes), use.ave.cell)
sqw <- rep(c(1, 1e-3), c(ncells*length(sizes), ncells))
fit <- limma::lmFit(design=as.matrix(collected[[1]])/sqw, object=collected[[2]]/sqw, weight=sqw^2) # Checking to limma's value.
expect_equal(unname(attributes(outx)$standard.error), unname(fit$sigma * fit$coefficients[1,]/as.numeric(outx)))

## obs.sf <- rnorm(200, mean=5) # True size factors for all cells (not quite realistic, as we shouldn't get negative values)..
## pooled <- as.matrix(collected[[1]]) %*% obs.sf
## pooled <- pooled + rnorm(length(pooled), sd=0.5) # estimation error in each pool (not quite realistic, as error scales with the mean; but oh well).
## fit <- limma::lmFit(design=as.matrix(collected[[1]]), object=t(pooled))
## fit$sigma # Should be the estimation error, equal to 0.5.
## head(as.numeric(fit$sigma  * fit$stdev.unscaled)) # NOT the estimation error as it is huge.

# Trying it out on a SCESet object.

count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
out <- computeSumFactors(X)
expect_equal(unname(sizeFactors(out)), computeSumFactors(dummy))

# Throwing in some silly inputs.

expect_error(computeSumFactors(dummy[,0,drop=FALSE]), "not enough cells in each cluster")
expect_error(computeSumFactors(dummy[0,,drop=FALSE]), "cells should have non-zero library sizes")
expect_warning(computeSumFactors(dummy[,1:100]), "number of cells in each cluster should be at least twice")
expect_error(computeSumFactors(dummy, sizes=c(10, 10, 20)), "'sizes' is not unique")
expect_error(computeSumFactors(dummy, clusters=integer(0)), "'x' ncols is not equal to 'clusters' length")

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

# Checking the matrix construction.

cur.exprs <- matrix(1, nrow=ngenes, ncol=ncells)
subsphere <- sample(ncells)
sphere <- c(subsphere, subsphere)
sizes <- c(20, 51)
use.ave.cell <- rep(1, ngenes)
out <- scran:::.create_linear_system(cur.exprs, sphere, as.integer(sizes), use.ave.cell)

out$design <- as.matrix(out$design)
size <- sizes[1]
for (i in seq_len(ncells)) { 
    used <- sphere[i:(i+size-1)]
    expect_equal(out[[1]][i, used], rep(1, size))
    expect_equal(out[[1]][i, -used], rep(0, ncells-size))
}
size <- sizes[2]
for (i in seq_len(ncells)) { 
    used <- sphere[i:(i+size-1)]
    expect_equal(out[[1]][i + ncells, used], rep(1, size))
    expect_equal(out[[1]][i + ncells, -used], rep(0, ncells-size))
}
for (i in seq_len(ncells)) { 
    used <- sphere[i]
    expect_equal(out[[1]][i + ncells*2, used], sqrt(0.000001))
    expect_equal(out[[1]][i + ncells*2, -used], rep(0, ncells-1L))
}

expect_identical(nrow(out$design), as.integer(ncells*(length(sizes)+1L)))
expect_identical(ncol(out$design), as.integer(ncells))
expect_equal(out[[2]], rep(c(sizes, sqrt(0.000001)), each=ncells))

# Checking some other internals.

x <- matrix(rpois(20000, lambda=10), nrow=200, ncol=100)
subset.row <- sample(nrow(x), 50)
subset.col <- sample(ncol(x), 50)
cur.out <- .Call(scran:::cxx_subset_and_divide, x, subset.row-1L, subset.col-1L)
expect_equal(cur.out[[1]], colSums(x[subset.row,subset.col]))
expect_equal(cur.out[[2]], t(t(x[subset.row,subset.col])/cur.out[[1]]))

####################################################################################################

# Checking out what happens with clustering.

set.seed(20001)
ncells <- 700
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- seq_len(ngenes)/100
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0

emp.clusters <- quickCluster(dummy)
expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
shuffled <- c(1:50, 301:350, 601:650)
expect_identical(quickCluster(dummy, subset.row=shuffled), emp.clusters)

# Checking out deeper internals.

set.seed(200011)
mat <- matrix(rpois(10000, lambda=5), nrow=20)

subset.row <- seq_len(nrow(mat))
distM <- .Call(scran:::cxx_compute_cordist, mat, subset.row - 1L)
refM <- sqrt(0.5*(1 - cor(mat, method="spearman")))
expect_equal(distM, refM)

subset.row <- 15:1 # With subsetting
distM <- .Call(scran:::cxx_compute_cordist, mat, subset.row - 1L)
refM <- sqrt(0.5*(1 - cor(mat[subset.row,], method="spearman")))
expect_equal(distM, refM)

suppressWarnings(expect_false(identical(quickCluster(mat), quickCluster(mat[subset.row,])))) # Checking that subsetting gets different results.
suppressWarnings(expect_identical(quickCluster(mat, subset.row=subset.row), quickCluster(mat[subset.row,]))) # Checking that subset.row works.

# Checking out what happens with silly inputs.

expect_error(quickCluster(dummy[0,]), "need at least 2 observations to compute correlations")
expect_error(quickCluster(dummy[,0]), "fewer cells than the mininimum cluster size")

leftovers <- 100
expect_warning(forced <- quickCluster(dummy[,c(which(known.clusters==1), which(known.clusters==2), which(known.clusters==3)[1:leftovers])]), 
               sprintf("%i cells were not assigned to any cluster", leftovers))
expect_identical(as.character(tail(forced, leftovers)), rep("0", leftovers))

# Seeing how it interacts with the normalization method.

out <- computeSumFactors(dummy, cluster=known.clusters)
expect_equal(out, count.sizes/mean(count.sizes)) # Even though there is a majority of DE, each pair of clusters is still okay.

out1 <- computeSumFactors(dummy, cluster=known.clusters, ref=1)
expect_equal(out, out1)
out2 <- computeSumFactors(dummy, cluster=known.clusters, ref=2)
expect_equal(out, out2)
out3 <- computeSumFactors(dummy, cluster=known.clusters, ref=3)
expect_equal(out, out3)

expect_error(computeSumFactors(dummy, cluster=known.clusters, ref=0), "'ref.clust' value not in 'clusters'")

# Trying it out on a SCESet object.

set.seed(20002)
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- sample(seq_len(ngenes)/100)
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
emp.clusters <- quickCluster(X)
expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)

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
isSpike(X) <- "MySpike"
out <- computeSpikeFactors(X)
ref <- colSums(dummy[is.spike,])
expect_equal(unname(sizeFactors(out)), ref/mean(ref))
expect_equal(sizeFactors(out), sizeFactors(out, type="MySpike"))

# Checking out what happens when you have multiple spike-ins supplied.
X2 <- newSCESet(countData=data.frame(dummy))
X2 <- calculateQCMetrics(X2, list(MySpike=is.spike, SecondSpike=is.spike))
isSpike(X2) <- c("MySpike", "SecondSpike")
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
sf <- ref/mean(ref)
sizeFactors(X) <- sf
out <- normalize(X)
expect_equivalent(exprs(out), edgeR::cpm.default(dummy, lib.size=sf*1e6, prior.count=1, log=TRUE))
out <- normalize(X, logExprsOffset=3)
expect_equivalent(exprs(out), edgeR::cpm.default(dummy, lib.size=sf*1e6, prior.count=3, log=TRUE))

sf <- runif(ncells, 10, 20)
sf <- ref/mean(ref)
sizeFactors(X) <- sf
out <- normalize(X)
expect_equivalent(exprs(out), edgeR::cpm.default(dummy, lib.size=sf*1e6, prior.count=1, log=TRUE))

chosen <- rbinom(ngenes, 1, 0.7)==0L
X <- calculateQCMetrics(X, feature_controls=list(whee=chosen))
X3 <- normalize(X)
expect_equivalent(exprs(out), exprs(X3))

sizeFactors(X, type="whee") <- colSums(counts(X)[chosen,])
X4 <- normalize(X)
expect_equivalent(exprs(out)[!chosen,], exprs(X4)[!chosen,])
ref <- sizeFactors(X, type="whee")
sf <- ref/mean(ref)
expect_equivalent(exprs(X4)[chosen,], edgeR::cpm.default(dummy[chosen,], lib.size=sf*1e6, prior.count=1, log=TRUE))

# Checking out silly inputs.

expect_equal(unname(dim(normalize(X[,0,drop=FALSE]))), c(ngenes, 0L))
expect_equal(unname(dim(normalize(X[0,,drop=FALSE]))), c(0L, ncells)) 


