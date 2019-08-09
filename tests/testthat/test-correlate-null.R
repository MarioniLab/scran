# This checks the correlateNull function.
# require(scran); require(testthat); source("setup.R"); source("test-correlate-null.R")

set.seed(20000)
test_that("null distribution of correlations looks okay", {
    for (ncells in c(20, 50, 100)) {
        out <- correlateNull(ncells, iters=1e3)
        expect_equal(length(out), 1e3)
        expect_lte(max(out), 1)
        expect_gte(min(out), -1)

        random1 <- matrix(rnorm(ncells*1000), ncol=1000)
        random2 <- matrix(rnorm(ncells*1000), ncol=1000)
        ref <- as.vector(cor(random1, random2, method="spearman"))

        expect_true(abs(mean(out) - mean(ref)) < 0.01)
        expect_true(abs(var(out) - var(ref)) < 0.01)
    }

    # Responds to a seed.
    set.seed(0)
    out1 <- correlateNull(121, iters=1e3)
    set.seed(0)
    out2 <- correlateNull(121, iters=1e3)
    expect_identical(out1, out2)
})

set.seed(20001)
test_that("correlateNull works with a blocking factor", {
    grouping <- sample(LETTERS[1:5], 121, replace=TRUE)

    set.seed(100)
    out1 <- out2 <- 0
    for (gl in table(grouping)) {
        X <- correlateNull(gl, 1e3)
        out1 <- out1 + X * gl
        out2 <- out2 + X
    }
    out1 <- out1/length(grouping)
    out2 <- out2/length(unique(grouping))

    set.seed(100)
    out1x <- correlateNull(block=grouping, equiweight=FALSE, iters=1e3)
    set.seed(100)
    out2x <- correlateNull(block=grouping, iters=1e3)

    expect_equal(out1, out1x)
    expect_equal(out2, out2x)

    # Ignores blank blocks.
    set.seed(200)
    ref <- correlateNull(block=grouping, iters=1e3)
    set.seed(200)
    out <- correlateNull(block=c("Z", grouping), iters=1e3)
    expect_equal(out, ref)
    set.seed(200)
    out <- correlateNull(block=c("Z", "Z", grouping), iters=1e3)
    expect_equal(out, ref)
})

set.seed(20002)
test_that("correlateNull works with a design matrix", {
    for (design in list(
        oneway=model.matrix(~factor(rep(c(1,2), each=10))),
        covariate=model.matrix(~seq_len(10))
        ))
    {
        out <- correlateNull(design=design, iters=1e5)

        ncells <- nrow(design)
        random1 <- lm.fit(y=matrix(rnorm(ncells*1000), ncol=1000), x=design)$residuals
        random2 <- lm.fit(y=matrix(rnorm(ncells*1000), ncol=1000), x=design)$residuals
        ref <- as.vector(cor(random1, random2, method="spearman"))

        expect_true(abs(mean(out) - mean(ref)) < 0.01)
        expect_true(abs(var(out) - var(ref)) < 0.01)
    }
})

test_that("correlateNull is unaffected by the number of cores", {
    set.seed(200)
    ref <- correlateNull(12, iters=1e3)

    BPPARAM <- safeBPParam(3) # Before set.seed, as safeBPParam resets the seed.
    set.seed(200)
    out <- correlateNull(12, iters=1e3, BPPARAM=BPPARAM)
    expect_identical(ref, out)

    grouping <- sample(rep(LETTERS[1:5], each=5))
    set.seed(300)
    ref <- correlateNull(block=grouping, iters=1e3)
    set.seed(300)
    out <- correlateNull(block=grouping, iters=1e3, BPPARAM=BPPARAM)
    expect_identical(ref, out)

    design <- cbind(1, runif(30))
    set.seed(400)
    ref <- correlateNull(design=design, iters=1e3)
    set.seed(400)
    out <- correlateNull(design=design, iters=1e3, BPPARAM=BPPARAM)
    expect_identical(ref, out)

    # Does not fail with low numbers of iterations.
    expect_error(correlateNull(12, iters=2, BPPARAM=BPPARAM), NA)
})

test_that("correlateNull works correctly on silly inputs", {
    expect_identical(correlateNull(ncells=100, iters=0), numeric(0))
    expect_identical(correlateNull(ncells=0, iters=1), NaN)
    expect_identical(correlateNull(design=diag(10), iters=1), NaN)

    expect_error(correlateNull(200, block=rep(1, 20)), "cannot specify")
    expect_error(correlateNull(200, design=cbind(rep(1, 20))), "cannot specify")
})
