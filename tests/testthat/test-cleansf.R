# This tests the cleanSizeFactors() function.
# library(scran); library(testthat); source("test-cleansf.R")

set.seed(1000)
test_that("cleanSizeFactors works correctly", {
    N <- 100
    sf <- runif(N)
    num <- jitter(sf/(sf+1)) # nls() fails on zero-residual data.

    # Does nothing; just returns the size factors directly.
    out <- cleanSizeFactors(sf, num)
    expect_equal(sf, out)
    
    # Fills in the expected values. 
    expect_warning(out <- cleanSizeFactors(c(-sf, sf), c(num, num)))
    expect_identical(sf, tail(out, N))
    expect_equal(sf, head(out, N), tol=1e-3)

    # Handles cases where the number of genes is greater than possible.
    expect_warning(out <- cleanSizeFactors(c(-1, sf), c(2, num)))
    expect_equal(out[1], max(sf))

    # Avoids warnings when there is actually a decent amount of noise.
    jnum <- num * 2^rnorm(N, sd=0.1) 
    expect_warning(out <- cleanSizeFactors(c(-1, sf), c(0.5, jnum)), NA)
    expect_identical(sf, tail(out, -1))
})

set.seed(1001)
test_that("robustness iterations make a difference in cleanSizeFactors", {
    N <- 100
    sf <- runif(N)
    num <- jitter(sf/(sf+1))

    sf2 <- c(sf, 0.1, -sf) # adding an outlier.
    num2 <- c(num, 1, num)

    expect_warning(out <- cleanSizeFactors(sf2, num2))
    expect_identical(sf, head(out, N))
    expect_equal(sf, tail(out, N), tol=1e-3)

    # Without robustness, it does not return accurate values.
    out <- cleanSizeFactors(sf2, num2, iterations=0)
    expect_identical(sf, head(out, N))
    expect_false(all(abs(tail(out, N) - sf) < 1e-3))
})

set.seed(1002)
test_that("cleanSizeFactors behaves in a live example", {
    counts <- matrix(rpois(20000, lambda=runif(100)), ncol=100, byrow=TRUE)

    # Adding negative values:
    num <- c(100, colSums(counts>0))
    sf <- c(-1, colSums(counts))
    out <- cleanSizeFactors(sf, num)
    
    expect_identical(tail(sf, -1), tail(out, -1))
    expect_true(out[1] > mean(sf[num < 100]))
    expect_true(out[1] < mean(sf[num > 100]))
})

test_that("cleanSizeFactors avoids silly inputs", {
    expect_identical(cleanSizeFactors(numeric(0), integer(0)), numeric(0))        
    expect_error(cleanSizeFactors(-1, integer(0)), "same length")        
    expect_error(cleanSizeFactors(-(1:10), 1:10), "need at least")        
})
