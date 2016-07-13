# This tests the various SCESet methods in scran.

require(scran); require(testthat);

set.seed(30000)
ncells <- 200
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- newSCESet(countData=data.frame(dummy))
is.spike <- rbinom(ngenes, 1, 0.5)==0L
X <- calculateQCMetrics(X, feature_controls=list(whee=is.spike, empty=logical(ngenes), all=!logical(ngenes)))
isSpike(X) <- "whee"
expect_identical(isSpike(X), is.spike)
expect_identical(isSpike(X, type="whee"), is.spike)

sf <- runif(ncells, 0.5, 1.5)
sizeFactors(X) <- sf
expect_identical(sf, unname(sizeFactors(X)))
expect_identical(colnames(X), names(sizeFactors(X)))

expect_identical(spikes(X), counts(X)[isSpike(X),,drop=FALSE])
X <- normalize(X)
expect_identical(spikes(X, "exprs"), exprs(X)[isSpike(X),,drop=FALSE])

isSpike(X) <- c("whee", "empty") # Testing multiple specifications.
expect_identical(isSpike(X), is.spike)
isSpike(X) <- c("whee", "all") 
expect_identical(isSpike(X), !logical(ngenes))
isSpike(X) <- c("whee", "empty", "all") 
expect_identical(isSpike(X), !logical(ngenes))

isSpike(X) <- NULL
expect_warning(isSpike(X), "'isSpike' is not set, returning NULL")

# Checking silly inputs

isSpike(X) <- "whee"
sizeFactors(X) <- sf
expect_error(isSpike(X) <- "aaron", "'aaron' is not specified as a spike-in control")
expect_error(sizeFactors(X) <- "whee", "unable to find an inherited method")
expect_identical(isSpike(X[0,]), logical(0))
expect_identical(unname(sizeFactors(X[,0])), numeric(0))

expect_identical(spikes(X[0,]), exprs(X)[0,])
expect_identical(spikes(X[,0]), exprs(X)[isSpike(X),0])
isSpike(X) <- "empty"
expect_identical(spikes(X), exprs(X)[0,])

isSpike(X) <- NULL
expect_warning(out <- isSpike(X), "'isSpike' is not set, returning NULL")
expect_identical(out, NULL)

sizeFactors(X) <- NULL
expect_warning(out <- sizeFactors(X), "not been set")
expect_identical(out, NULL)
