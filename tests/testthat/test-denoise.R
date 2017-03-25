# This checks the denoisePCA function.

# require(scran); require(testthat); source("test-denoise.R")

# Mocking up some data with subpopulations of cells.

set.seed(1000)
ngenes <- 1000
npops <- 5
ncells <- 100
means <- 2^runif(ngenes, 6, 10)
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

# Checking that the output is the same.

pcs <- denoisePCA(lcounts, technical=fit$trend)
pcs2 <- denoisePCA(lcounts, technical=setNames(dec$tech, rownames(dec)))
expect_equal(pcs, pcs2)

pcs3 <- denoisePCA(lcounts, technical=fit$trend, design=cbind(rep(1, ncells)))
expect_equal(pcs, pcs3)

not.spike <- setdiff(seq_len(ngenes), is.spike)
pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
pcs2 <- denoisePCA(lcounts[not.spike,], technical=fit$trend)
expect_equal(pcs, pcs2)

# Checking for sensible handling of design matrices.

design <- model.matrix(~factor(in.pop))
dfit <- trendVar(lcounts, subset.row=is.spike, design=design)
pcs <- denoisePCA(lcounts, design=design, technical=dfit$trend)

alt <- lm.fit(y=t(lcounts), x=design)
true.var <- colMeans(alt$effects[-seq_len(alt$rank),]^2)
obs.var <- apply(alt$residuals, 2, var)
new.x <- alt$residuals * sqrt(true.var/obs.var)
alt.pc <- prcomp(new.x)
expect_equal(alt.pc$x[,seq_len(ncol(pcs)),drop=FALSE], pcs)

# Checking invalid specifications.

expect_error(denoisePCA(lcounts, technical=c(Whee=1)), "missing gene names in 'technical'")
unnamed.lcounts <- lcounts
rownames(unnamed.lcounts) <- NULL
expect_error(denoisePCA(unnamed.lcounts, technical=c(Whee=1)), "rows of 'x' should be named with gene names")

# Checking for proper behaviour with SCESet.

X <- newSCESet(exprsData=lcounts, logExprsOffset=1, lowerDetectionLimit=0)
X2 <- denoisePCA(X, technical=fit$trend)
pcx <- reducedDimension(X2)
rownames(pcx) <- NULL
pcs <- denoisePCA(lcounts, technical=fit$trend)
expect_equal(pcx, pcs)

X <- calculateQCMetrics(X, feature_controls=list(Spike=is.spike))
setSpike(X) <- "Spike"
X2 <- denoisePCA(X, technical=fit$trend, get.spikes=TRUE)
pcx <- reducedDimension(X2)
rownames(pcx) <- NULL
expect_equal(pcx, pcs)

X2 <- denoisePCA(X, technical=fit$trend)
pcx <- reducedDimension(X2)
rownames(pcx) <- NULL
pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=not.spike)
expect_equal(pcx, pcs)


