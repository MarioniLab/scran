# This tests the various trendVar() options.
# require(scran); require(testthat); source("test-trend-var.R")

set.seed(20000)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

library(scater)
out <- normalizeCounts(dummy)

library(DelayedMatrixStats)
means <- rowMeans(out)
vars <- rowVars(out)

test_that("trendVar works on a basic scenario", {
    out <- fitTrendVar(means, vars)
    expect_true(out$std.dev > 0)
    expect_is(out$trend, "function")

    # Hard to test it without copying the code, so I'll just check the limits.
    expect_equal(out$trend(0), 0)
    expect_equal(out$trend(1:10), sapply(1:10, out$trend)) # Checking we get consistent results with returned function.
    expect_equal(out$trend(100:1/20), sapply(100:1/20, out$trend))  # More checking, reversed order.

    # Parametric left edge limits work (almost) correctly (below 0.1).
    expect_equal(out$trend(0.01)*2, out$trend(0.02), tol=1e-4)
    expect_equal(out$trend(0.01)*5, out$trend(0.05), tol=1e-4)

    # Parametric right limits do not use rule=2.
    mx <- max(means)
    expect_false(isTRUE(all.equal(out$trend(mx), out$trend(mx+1))))
})

test_that("trendVar works when parametric mode is turned off", {
    mx <- max(means)
    out <- fitTrendVar(means, vars, parametric=FALSE)
    expect_equal(out$trend(0), 0)
    expect_equal(out$trend(1:10), sapply(1:10, out$trend)) 
    expect_equal(out$trend(100:1/20), sapply(100:1/20, out$trend)) 

    # Parametric left edge limits work correctly (below 0.1).
    expect_equal(out$trend(0.01)*2, out$trend(0.02))
    expect_equal(out$trend(0.01)*5, out$trend(0.05))
    expect_equal(out$trend(0.01)*10, out$trend(0.1))

    # Right edge goes to rule=2.
    expect_equal(out$trend(mx), out$trend(mx+1))
})

test_that("trendVar handles nls errors gracefully",  {
    X <- runif(100)
    Y <- runif(100)
    expect_warning(out <- fitTrendVar(X, Y), "defaulting to loess-only")
    expect_warning(ref <- fitTrendVar(X, Y, parametric=FALSE), NA)
    expect_identical(out$std.dev, ref$std.dev)

    expect_error(fitTrendVar(runif(2), runif(2)), "need at least 4") 
    expect_error(fitTrendVar(runif(1), runif(1)), "need at least 2") 
})