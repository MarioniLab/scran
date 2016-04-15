# This tests whether the shuffling procedure is doing its job.

set.seed(0)

N <- 1000L
rankings <- as.double(seq_len(3))
collected <- list()
for (it in 1:500) {
    my.shuffle <- .Call(scran:::cxx_auto_shuffle, rankings, N)
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
expect_true(all(abs(colMeans(out) - N/length(rankings)) < 2))

