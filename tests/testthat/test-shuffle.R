# This tests whether the shuffling procedure is doing its job.
# library(testthat); source("setup.R"); source("test-shuffle.R")

test_that("shuffling is behaving correctly", {
    set.seed(0)
    
    N <- 1000L
    rankings <- as.double(seq_len(3))
    collected <- list()
    for (it in 1:500) {
        my.shuffle <- scrambler(rankings, N, reset=FALSE)
        expect_identical(nrow(my.shuffle), length(rankings))
        expect_identical(ncol(my.shuffle), N)
    
        is.1 <- rowSums(my.shuffle==1)
        is.2 <- rowSums(my.shuffle==2)
        is.3 <- rowSums(my.shuffle==3)
        expect_identical(as.integer(sum(is.1)), N)
        expect_identical(as.integer(sum(is.2)), N)
        expect_identical(as.integer(sum(is.3)), N)
        expect_identical(as.integer(is.1 + is.2 + is.3), rep(N, 3))
    
        collected[[it]] <- c(is.1, is.2, is.3)
    }
    
    out <- do.call(rbind, collected)
    expect_true(all(abs(colMeans(out) - N/length(rankings)) < 2)) # Should be very close to the expectation.
})

test_that("shuffling behaves with reset=TRUE", {
    set.seed(0)
    
    N <- 1000L
    rankings <- as.double(seq_len(3))
    collected <- list()
    for (it in 1:500) {
        my.shuffle <- scrambler(rankings, N, reset=TRUE)
        expect_identical(nrow(my.shuffle), length(rankings))
        expect_identical(ncol(my.shuffle), N)
    
        is.1 <- rowSums(my.shuffle==1)
        is.2 <- rowSums(my.shuffle==2)
        is.3 <- rowSums(my.shuffle==3)
        expect_identical(as.integer(sum(is.1)), N)
        expect_identical(as.integer(sum(is.2)), N)
        expect_identical(as.integer(sum(is.3)), N)
        expect_identical(as.integer(is.1 + is.2 + is.3), rep(N, 3))
    
        collected[[it]] <- c(is.1, is.2, is.3)
    }
    
    out <- do.call(rbind, collected)
    expect_true(all(abs(colMeans(out) - N/length(rankings)) < 2)) # Should be very close to the expectation.
})

test_that("shuffling is responding to the seed", {
    set.seed(0)
    for (seed in 10^0:5) {
        blah <- runif(100)
        N <- 20

        set.seed(seed)
        out1 <- scrambler(blah, N)
        out2 <- scrambler(blah, N)
        expect_false(all(out1==out2)) # Should be different.

        set.seed(seed)
        out3 <- scrambler(blah, N)
        expect_identical(out1, out3) # Should be the same.
    }
})
