# This tests the makeTechTrend function.
# require(scran); require(testthat); source("test-techtrend.R")

REFFUN <- function(means, size.factors, tol=1e-6, dispersion=0, pseudo.count=1) {
    # Defining which distribution to use.
    if (dispersion==0) {
        qfun <- function(..., mean) { qpois(..., lambda=mean) }
        dfun <- function(..., mean) { dpois(..., lambda=mean) }
    } else {
        qfun <- function(..., mean) { qnbinom(..., mu=mean, size=1/dispersion) }
        dfun <- function(..., mean) { dnbinom(..., mu=mean, size=1/dispersion) }
    }
    collected.means <- collected.vars <- numeric(length(means))

    for (i in seq_along(means)) {
        m <- means[i]*size.factors
        lower <- qfun(tol, mean=m, lower.tail=TRUE)
        upper <- qfun(tol, mean=m, lower.tail=FALSE)

        # Creating a function to compute the relevant statistics.
        .getValues <- function(j) {
            ranged <- lower[j]:upper[j]
            p <- dfun(ranged, mean=m[j])
            lvals <- log2(ranged/size.factors[j] + pseudo.count)
            return(list(p=p, lvals=lvals))
        }

        # Computing the mean.
        cur.means <- numeric(length(size.factors))
        for (j in seq_along(size.factors)) { 
            out <- .getValues(j)
            cur.means[j] <- sum(out$lvals * out$p) / sum(out$p)
        }
        final.mean <- mean(cur.means)
        collected.means[i] <- final.mean
       
        # Computing the variance. Done separately to avoid
        # storing 'p' and 'lvals' in memory, but as a result
        # we need to compute these values twice.
        cur.vars <- numeric(length(size.factors))
        for (j in seq_along(size.factors)) { 
            out <- .getValues(j)
            cur.vars[j] <- sum((out$lvals - final.mean)^2 * out$p) / sum(out$p)
        }
        collected.vars[i] <- mean(cur.vars)
    }

    return(list(mean=collected.means, var=collected.vars))
}

set.seed(20004)
test_that("makeTechTrend compares correctly to a reference", {
    # Checking the C++ function against a reference R implementation.
    sf <- runif(100)
    means <- sort(runif(20, 0, 50))
    ref <- REFFUN(means, sf, tol=1e-6, dispersion=0, pseudo.count=1)
    out <- scran:::calc_log_count_stats(means, sf, 1e-6, 0, 1)
    expect_equal(ref[[1]], out[[1]])
    expect_equal(ref[[2]], out[[2]])

    ref <- REFFUN(means, sf, tol=1e-6, dispersion=0.1, pseudo.count=1)
    out <- scran:::calc_log_count_stats(means, sf, 1e-6, 0.1, 1)
    expect_equal(ref[[1]], out[[1]])
    expect_equal(ref[[2]], out[[2]])
})

test_that("makeTechTrend values are close to simulated means/variances", {
    logrpois <- function(m, s, p) {
        log2(rpois(1e6, lambda=m*s)/s + p)
    }

    out <- makeTechTrend(c(1, 5), pseudo.count=2, size.factors=c(0.5, 1.5))
    log.values <- c(logrpois(1, 0.5, 2), logrpois(1, 1.5, 2)) # at 1
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)
    log.values <- c(logrpois(5, 0.5, 2), logrpois(5, 1.5, 2)) # at 5
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)

    # With a NB distribution.
    out <- makeTechTrend(c(1, 5), dispersion=0.1)
    log.values <- log2(rnbinom(1e6, mu=1, size=10) + 1)
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)
    log.values <- log2(rnbinom(1e6, mu=5, size=10) + 1)
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)
})

test_that("makeTechTrend performs approximations correctly", {
    sf <- runif(200)
    sf <- sf/mean(sf)

    # Checking that the calculations are exactly correct with perfect interpolation.
    out <- makeTechTrend(1:100/10, pseudo.count=1, size.factors=sf)
    out2 <- makeTechTrend(1:100/10, pseudo.count=1, size.factors=sf, approx.npts=200)
    expect_equal(environment(out)$z$x, environment(out2)$z$x)
    expect_equal(environment(out)$z$y, environment(out2)$z$y)

    # Checking that the calculations are close enough with non-trivial interpolation.
    out <- makeTechTrend(1:100/10, pseudo.count=1, size.factors=sf)
    out2 <- makeTechTrend(1:100/10, pseudo.count=1, size.factors=sf, approx.npts=100)
    expect_equal(environment(out)$z$x, environment(out2)$z$x, 1e-6)
    expect_equal(environment(out)$z$y, environment(out2)$z$y, 1e-6)

    # Throws errors.
    expect_error(makeTechTrend(1:100, pseudo.count=1, size.factors=sf, approx.npts=1), "at least 2")
})

test_that("makeTechTrend handles other options properly", {
    # Handles parallelization properly.
    sf <- c(0.1, 0.5, 1, 1.5, 1.9)
    means <- 0:20/5
    out1 <- makeTechTrend(means, pseudo.count=1, size.factors=sf)
    out2 <- makeTechTrend(means, pseudo.count=1, size.factors=sf, BPPARAM=safeBPParam(2))
    expect_equal(out1(means), out2(means))

    # Returns a zero when a zero is provided.
    expect_equal(out1(0), 0)

    # Happy with non-centred size factors.
    expect_error(makeTechTrend(0:5, size.factors=1:5), NA)
})        

test_that("makeTechTrend handles SCE inputs correctly", {
    X <- SingleCellExperiment(list(counts=matrix(1:5, ncol=2, nrow=5)))
    expect_error(makeTechTrend(x=X), "log.exprs.offset")

    sizeFactors(X) <- c(0.9, 1.1)
    suppressWarnings(X <- normalize(X))
    out <- makeTechTrend(x=X)
    ref <- makeTechTrend(2^seq(0, max(rowMeans(exprs(X))), length.out=100)-1,
                         size.factors=sizeFactors(X))
    expect_equal(out(0:10/2), ref(0:10/2))

    # Automatically uses the library sizes as size factors.
    X <- SingleCellExperiment(list(counts=matrix(1:10, ncol=2, nrow=5)))
    suppressWarnings(X <- normalize(X))
    out <- makeTechTrend(x=X)
    libsizes <- colSums(counts(X))
    ref <- makeTechTrend(2^seq(0, max(rowMeans(exprs(X))), length.out=100)-1,
                         size.factors=libsizes/mean(libsizes))
    expect_equal(out(0:10/2), ref(0:10/2))

    # Chucks an warning when there are multiple features with non-centred size factors.
    isSpike(X, "ERCC") <- 1
    sizeFactors(X) <- runif(ncol(X))
    sizeFactors(X, "ERCC") <- runif(ncol(X))
    expect_warning(makeTechTrend(x=X), "centred") 
})

