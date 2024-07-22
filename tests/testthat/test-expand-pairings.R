# This tests the expand_pairings utility in scran.
# library(testthat); library(scran); source("test-expand-pairings.R")

test_that("expand_pairings works as expected for NULL", {
    out <- scran:::.expand_pairings(NULL, universe=1:5)
    expect_identical(out$mode, "single pool")

    grid <- expand.grid(1:5, 1:5)
    grid <- grid[grid[,1] != grid[,2],]
    expect_identical(out$id1, grid[,1])
    expect_identical(out$id2, grid[,2])

    # Even when empty.
    out <- scran:::.expand_pairings(NULL, universe=integer(0))
    expect_identical(out$id1, integer(0))
    expect_identical(out$id2, integer(0))
})

test_that("expand_pairings works as expected for vectors", {
    out <- scran:::.expand_pairings(2:4, universe=1:5)
    expect_identical(out$mode, "single pool")

    grid <- expand.grid(2:4, 2:4)
    grid <- grid[grid[,1] != grid[,2],]
    expect_identical(out$id1, grid[,1])
    expect_identical(out$id2, grid[,2])

    out2 <- scran:::.expand_pairings(c("C", "D", "E"), universe=LETTERS[2:5])
    expect_identical(out, out2)

    ref2 <- scran:::.expand_pairings(1:3, universe=1:5)
    out2 <- scran:::.expand_pairings(5:3, universe=5:1)
    expect_identical(ref2, out2)

    # Works for elements out of range.
    out3 <- scran:::.expand_pairings(c("C", "D", "E"), universe=LETTERS[1:4])
    out4 <- scran:::.expand_pairings(3:4, universe=1:4)
    expect_identical(out3, out4)

    # Even when empty.
    out <- scran:::.expand_pairings(integer(0), universe=1:5)
    expect_identical(out$id1, integer(0))
    expect_identical(out$id2, integer(0))
})

test_that("expand_pairings works as expected for lists", {
    out <- scran:::.expand_pairings(list(c(1,5), c(2:4)), universe=1:5)
    expect_identical(out$mode, "double pool")

    grid <- expand.grid(c(1L,5L), 2:4)
    expect_identical(out$id1, grid[,1])
    expect_identical(out$id2, grid[,2])

    ref2 <- scran:::.expand_pairings(list(1:3, 4:5), universe=1:5)
    out2 <- scran:::.expand_pairings(list(c(5,3,1),c(4,2)), universe=c(5,3,1,4,2))
    expect_identical(ref2, out2)

    # Works for character vectors.
    out2 <- scran:::.expand_pairings(list(c("B", "F"), c("C", "D", "E")), universe=LETTERS[2:6])
    expect_identical(out, out2)

    # Works for NULL.
    null <- scran:::.expand_pairings(list(NULL, NULL), universe=LETTERS[2:6])
    ref <- scran:::.expand_pairings(NULL, universe=LETTERS[2:6])
    expect_identical(null[1:2], ref[1:2])

    # Works for elements out of range.
    out3 <- scran:::.expand_pairings(list(c(1,3,5), c(2, 4, 6)), universe=1:4)
    out4 <- scran:::.expand_pairings(list(c(1,3), c(2, 4)), universe=1:4)
    expect_identical(out3, out4)

    # Even when empty.
    out <- scran:::.expand_pairings(list(1:5, 1:5), universe=integer(0))
    expect_identical(out$id1, integer(0))
    expect_identical(out$id2, integer(0))
})

test_that("expand_pairings works as expected for matrices", {
    mat <- cbind(1:5, 2:6)
    out <- scran:::.expand_pairings(mat, universe=1:6)
    expect_identical(out$mode, "predefined pairs")
    expect_identical(out$id1, mat[,1])
    expect_identical(out$id2, mat[,2])

    # Filters out invalid pairs.
    ref <- scran:::.expand_pairings(mat[1:4,], universe=1:5)
    out <- scran:::.expand_pairings(mat, universe=1:5)
    expect_identical(ref, out)

    # Handles empty inputs.
    out <- scran:::.expand_pairings(mat[0,], universe=integer(0))
    expect_identical(out$id1, integer(0))
    expect_identical(out$id2, integer(0))
})
