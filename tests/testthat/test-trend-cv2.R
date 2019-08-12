# This tests the various trendCV2() options.
# require(scran); require(testthat); source("test-trend-cv2.R")

set.seed(20002)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

library(DelayedMatrixStats)
means <- rowMeans(dummy)
cv2s <- rowVars(dummy)/means^2

test_that("trendCV2 works on a basic scenario", {
    out <- fitTrendCV2(means, cv2s, ncells)
    expect_true(out$std.dev > 0)
    expect_is(out$trend, "function")

    # Hard to test it without copying the code, so I'll just check the limits.
    expect_equal(out$trend(0), Inf)
    expect_equal(out$trend(1:10), sapply(1:10, out$trend)) # Checking we get consistent results with returned function.
    expect_equal(out$trend(100:1/20), sapply(100:1/20, out$trend))  # More checking, reversed order.
})

test_that("trendCV2 prunes out capped CV2 values", {
    dummy[1:10,] <- 0
    dummy[1:10,1] <- 1:10*100
    means <- rowMeans(dummy)
    cv2s <- rowVars(dummy)/means^2

    out <- fitTrendCV2(means, cv2s, ncells)
    ref <- fitTrendCV2(means[-(1:10)], cv2s[-(1:10)], ncells)
   
    expect_identical(out$std.dev, ref$std.dev)
    expect_identical(out$trend(1:100), ref$trend(1:100))
})

test_that("trendCV2 handles nls errors gracefully",  {
    cv2.err <- jitter(10/1:10)
    m.err <- 1:10
    expect_error(fitTrendCV2(m.err, cv2.err, 10, simplified=FALSE), "singular gradient")
    expect_error(fit <- fitTrendCV2(m.err, cv2.err, 10), NA)
})
