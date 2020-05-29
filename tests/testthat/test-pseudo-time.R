# This tests the testPseudotime functionality.
# library(testthat); library(scran); source("test-pseudo-time.R")

test_that("basic tests work as expected", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- runif(100)

    out <- testPseudotime(y, u)
    expect_type(out$logFC, "double")

    out <- testPseudotime(y, u, get.spline.coef=TRUE)
    expect_type(out$spline1, "double")
})

test_that("multi-path tests work as expected", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- runif(100)
    
    p <- cbind(u, u)
    path1 <- rbinom(length(u), 1, 0.5)==0
    p[!path1,1] <- NA
    p[path1,2] <- NA

    out <- testPseudotime(y, p)
    expect_type(out$logFC.path1, "double")
    expect_type(out$logFC.path2, "double")

    sout <- testPseudotime(y, p, get.spline.coef=TRUE)
    expect_type(sout$spline1.path2, "double")

    # Removal of the cells in the shared path has no effect.
    p2 <- p
    p2[1:50,] <- runif(50)
    ref <- testPseudotime(y[,-(1:50)], p2[-(1:50),])
    out2 <- testPseudotime(y, p2)
    expect_identical(ref, out2)

    # Respects column names in 'p'.
    colnames(p) <- c("A", "B")

    out <- testPseudotime(y, p)
    expect_type(out$logFC.A, "double")
    expect_type(out$logFC.B, "double")

    out <- testPseudotime(y, p, get.spline.coef=TRUE)
    expect_type(out$spline1.B, "double")

    colnames(p) <- c("B", "B")
    expect_warning(testPseudotime(y, p), "duplicated")
})

test_that("handles large numbers of duplicated pseudotime values", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- c(numeric(50), runif(50))
    expect_error(out <- testPseudotime(y, u), NA)
    expect_error(out <- testPseudotime(y, u, df=51), "not enough unique")
})
