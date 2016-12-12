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
setSpike(X) <- "whee"
expect_identical(isSpike(X), is.spike)
expect_identical(whichSpike(X), "whee")
expect_identical(isSpike(X, type="whee"), is.spike)

sf <- runif(ncells, 0.5, 1.5)
sizeFactors(X) <- sf
expect_identical(sf, unname(sizeFactors(X)))
expect_identical(colnames(X), names(sizeFactors(X)))

expect_identical(spikes(X), counts(X)[isSpike(X),,drop=FALSE])
sizeFactors(X, "whee") <- sf # just to avoid the warning.
X <- normalize(X)
expect_identical(spikes(X, "exprs"), exprs(X)[isSpike(X),,drop=FALSE])

setSpike(X) <- c("whee", "empty") # Testing multiple specifications.
expect_identical(isSpike(X), is.spike)
expect_identical(whichSpike(X), c("whee", "empty"))
expect_warning(setSpike(X) <- c("whee", "all"), "overlapping spike-in sets detected")
expect_identical(isSpike(X), !logical(ngenes))
expect_warning(setSpike(X) <- c("whee", "empty", "all"), "overlapping spike-in sets detected")
expect_identical(isSpike(X), !logical(ngenes))

expect_identical(isSpike(X, type="whee"), is.spike)
expect_identical(isSpike(X, type=c("whee", "empty")), is.spike)
expect_identical(isSpike(X, type=c("whee", "all")), !logical(ngenes))
expect_identical(isSpike(X, type="empty"), logical(ngenes))

setSpike(X) <- NULL
expect_warning(isSpike(X), "no spike-ins specified, returning NULL")

# Checking silly inputs

setSpike(X) <- "whee"
sizeFactors(X) <- sf
expect_error(setSpike(X) <- "aaron", "'aaron' is not specified as a spike-in control")
expect_error(sizeFactors(X) <- "whee", "unable to find an inherited method")
expect_identical(isSpike(X[0,]), logical(0))
expect_identical(unname(sizeFactors(X[,0])), numeric(0))

expect_identical(spikes(X[0,]), exprs(X)[0,])
expect_identical(spikes(X[,0]), exprs(X)[isSpike(X),0])
setSpike(X) <- "empty"
expect_identical(spikes(X), exprs(X)[0,])

setSpike(X) <- NULL
expect_warning(out <- isSpike(X), "no spike-ins specified, returning NULL")
expect_identical(out, NULL)

sizeFactors(X) <- NULL
expect_warning(out <- sizeFactors(X), "not been set")
expect_identical(out, NULL)
