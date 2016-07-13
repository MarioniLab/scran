# Checks the technicalCV2 function.

require(scran); require(testthat)

set.seed(6000)
ngenes <- 10000
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), ncol=nsamples)

sf <- 2^rnorm(nsamples)
is.spike <- logical(ngenes)
is.spike[seq_len(500)] <- TRUE

# Testing the CV2 calculator. 
chosen <- which(is.spike)
stuff <- .Call(scran:::cxx_compute_CV2, counts, chosen - 1L, sf)
normed <- t(t(counts[chosen,])/sf)
expect_equal(stuff[[1]], rowMeans(normed))
expect_equal(stuff[[2]], apply(normed, 1, var))

chosen <- sample(ngenes, 1000)
stuff <- .Call(scran:::cxx_compute_CV2, counts, chosen - 1L, sf)
normed <- t(t(counts[chosen,])/sf)
expect_equal(stuff[[1]], rowMeans(normed))
expect_equal(stuff[[2]], apply(normed, 1, var))

# Comparing the SCESet and non-SCESet methods.
rownames(counts) <- paste0("X", seq_len(ngenes))
colnames(counts) <- paste0("Y", seq_len(nsamples))
X <- newSCESet(countData=counts)
X <- calculateQCMetrics(X, list(Spikes=is.spike))
isSpike(X) <- "Spikes"

sizeFactors(X) <- sf
sizeFactors(X, type="Spikes") <- 1

default <- technicalCV2(counts, is.spike, sf.cell=sf, sf.spike=rep(1, nsamples))
as.sceset <- technicalCV2(X, spike.type="Spikes")
expect_equal(default, as.sceset)

# Testing what happens when is.spike=NA.
all.used <- technicalCV2(counts, is.spike=NA)
counts2 <- rbind(counts, counts)
rownames(counts2) <- paste0("X", seq_len(nrow(counts2)))
is.spike2 <- rep(c(FALSE, TRUE), each=nrow(counts))
all.used2 <- technicalCV2(counts2, is.spike2)
expect_equal(all.used, all.used2[!is.spike2,])

all.used <- technicalCV2(counts(X), sf.cell=sizeFactors(X), sf.spike=sizeFactors(X), is.spike=NA)
all.used3 <- technicalCV2(X, spike.type=NA)
expect_equal(all.used, all.used3)

# Testing for silly inputs.
expect_error(technicalCV2(X, spike.type="whee"), "'arg' should be one of")
expect_error(technicalCV2(X[0,], spike.type="Spikes"), "need at least 2 spike-ins for trend fitting")
expect_error(technicalCV2(X[,0], spike.type="Spikes"), "need two or more cells to compute variances")

out <- technicalCV2(counts, is.spike=!logical(ngenes))
expect_true(all(is.na(out$p.value)))
expect_error(technicalCV2(counts, is.spike=logical(ngenes)), "need at least 2 spike-ins for trend fitting")

