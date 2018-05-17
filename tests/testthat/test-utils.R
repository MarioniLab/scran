# Checks various utility functions explicitly.
# require(scran); require(testthat); source("test-utils.R")

test_that("subset index conversion works correctly", {
    x <- matrix(0, 10, 20)

    # Integer
    expect_identical(scran:::.subset_to_index(1:5, x, byrow=FALSE), 1:5)
    expect_identical(scran:::.subset_to_index(1:5, x, byrow=TRUE), 1:5)

    expect_identical(scran:::.subset_to_index(11:20, x, byrow=FALSE), 11:20)
    expect_error(scran:::.subset_to_index(11:20, x, byrow=TRUE), "out of range")

    # Character
    y <- x
    rownames(y) <- LETTERS[1:10]
    expect_identical(scran:::.subset_to_index(LETTERS[1:5], y, byrow=TRUE), 1:5)
    expect_error(scran:::.subset_to_index(LETTERS[1:5], y, byrow=FALSE), "out of range")

    y <- x
    colnames(y) <- LETTERS[1:20]
    expect_identical(scran:::.subset_to_index(LETTERS[1:5], y, byrow=FALSE), 1:5)
    expect_error(scran:::.subset_to_index(LETTERS[1:5], y, byrow=TRUE), "out of range")

    # Logical
    set.seed(0)
    i <- rbinom(10, 1, 0.5)==0
    expect_identical(scran:::.subset_to_index(i, x, byrow=TRUE), which(i))
    i <- rbinom(20, 1, 0.5)==0
    expect_identical(scran:::.subset_to_index(i, x, byrow=FALSE), which(i))

    # Edge cases.
    expect_identical(scran:::.subset_to_index(numeric(0), x, byrow=FALSE), integer(0))
    expect_identical(scran:::.subset_to_index(character(0), x, byrow=FALSE), integer(0))
    expect_identical(scran:::.subset_to_index(logical(0), x, byrow=FALSE), integer(0))
})

test_that("worker assignment works correctly", {
    library(BiocParallel)
    wout <- scran:::.worker_assign(101, SerialParam())
    expect_identical(wout, list(seq_len(101)))

    wout <- scran:::.worker_assign(101, MulticoreParam(workers=2))
    expect_identical(length(wout), 2L) 
    expect_true(all(lengths(wout) >= floor(101/2)))
    expect_identical(unlist(wout), seq_len(101))

    wout <- scran:::.worker_assign(101, SnowParam(workers=3))
    expect_identical(length(wout), 3L) 
    expect_true(all(lengths(wout) >= floor(101/3)))
    expect_identical(unlist(wout), seq_len(101))

    # Checking vector splitting works.
    set.seed(10)
    thing <- runif(101)
    out <- scran:::.split_vector_by_workers(thing, wout)
    expect_identical(lengths(out), lengths(wout))
    expect_identical(unlist(out), thing)

    # Checking matrix splitting works.
    thing <- cbind(runif(101))
    out <- scran:::.split_matrix_by_workers(thing, wout, byrow=TRUE)
    expect_identical(vapply(out, FUN=nrow, FUN.VALUE=0L), lengths(wout))
    expect_identical(do.call(rbind, out), thing)

    thing <- rbind(runif(101))
    out <- scran:::.split_matrix_by_workers(thing, wout, byrow=FALSE)
    expect_identical(vapply(out, FUN=ncol, FUN.VALUE=0L), lengths(wout))
    expect_identical(do.call(cbind, out), thing)
})
