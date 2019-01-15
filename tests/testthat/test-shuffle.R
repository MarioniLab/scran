# This tests whether the shuffling procedure is doing its job.
# library(testthat); source("setup.R"); source("test-shuffle.R")

set.seed(100)
test_that("vector shuffling is behaving correctly", {
    N <- 1000L
    rankings <- as.double(seq_len(3))
    collected <- list()
    for (it in 1:500) {
        my.shuffle <- scramble_vector(rankings, N)
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

set.seed(101)
test_that("matrix shuffling behaves correctly", {
    N <- 1000L
    values <- scramble_vector(as.double(seq_len(3)), N)

    collected <- list()
    for (it in 1:500) {
        my.shuffle <- scramble_matrix(values)
        expect_identical(dim(values), dim(my.shuffle))

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
    expect_true(all(abs(colMeans(out) - N/nrow(values)) < 2)) # Should be very close to the expectation.
})

set.seed(102)
test_that("shuffling is responding to the seed", {
    for (seed in 10^0:5) {
        blah <- runif(100)
        N <- 20

        # Works for vector.
        set.seed(seed)
        out1 <- scramble_vector(blah, N)
        out2 <- scramble_vector(blah, N)
        expect_false(all(out1==out2)) # Should be different.

        sorted <- apply(out1, 2, sort)
        expect_identical(sorted, apply(out2, 2, sort))
        expect_identical(sorted, matrix(sort(blah), nrow=length(blah), ncol=N))

        set.seed(seed)
        out3 <- scramble_vector(blah, N)
        expect_identical(out1, out3) # Should be the same.

        # Works for matrix.
        whee <- matrix(rnorm(2000), ncol=20)

        set.seed(seed)
        out1 <- scramble_matrix(whee)
        out2 <- scramble_matrix(whee)
        expect_false(all(out1==out2)) # Should be different.

        sorted <- apply(out1, 2, sort)
        expect_identical(sorted, apply(out2, 2, sort))
        expect_identical(sorted, apply(whee, 2, sort))

        set.seed(seed)
        out3 <- scramble_matrix(whee)
        expect_identical(out1, out3) # Should be the same.
    }

    # Checking that the same seed with a different stream gives different results.
    blah <- rnorm(2000)
    expect_false(isTRUE(all.equal(scramble_vector(blah, 1, seed=0, stream=1), scramble_vector(blah, 1, seed=0, stream=2))))
    expect_true(identical(scramble_vector(blah, 1, seed=0, stream=1), scramble_vector(blah, 1, seed=0, stream=1))) # As a negative control
    expect_false(isTRUE(all.equal(scramble_vector(blah, 1, seed=0, stream=0), scramble_vector(blah, 1, seed=1, stream=0)))) # Positive control.

    # Checking that matrix shuffling is reproducible.
    whee <- matrix(rnorm(2000), ncol=20)
    seeds <- all_positive_integers(ncol(whee))
    ref <- scramble_matrix(whee, seed=seeds)
    sub.ref <- scramble_matrix(whee[,1:5], seed=seeds[1:5], stream=1:5)
    expect_identical(ref[,1:5], sub.ref)

    sub.ref <- scramble_matrix(whee[,11:5], seed=seeds[11:5], stream=11:5)
    expect_identical(ref[,11:5], sub.ref)
})
