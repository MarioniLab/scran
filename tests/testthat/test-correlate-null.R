# This checks the correlateNull function.
# require(scran); require(testthat); source("setup.R"); source("test-correlate-null.R")

refnull <- function(niters, ncells, resort=TRUE) {
    rankings <- as.double(seq_len(ncells))
    setup <- scran:::.setup_pcg_state(niters)
    shuffled <- scramble_matrix(matrix(rankings, nrow=ncells, ncol=niters), seed=setup$seeds[[1]], stream=setup$streams[[1]])
    out <- cor(shuffled, rankings, method="spearman")
    if (resort) { out <- sort(out) }
    out
}

set.seed(20000)
test_that("null distribution of correlations is correctly calculated", {
    set.seed(100)
    ref <- refnull(1e3, 121)
    set.seed(100)
    out <- correlateNull(121, iters=1e3)
    expect_equal(ref, as.double(out))
    
    set.seed(100)
    ref <- refnull(1e3, 12)
    set.seed(100)
    out <- correlateNull(12, iters=1e3)
    expect_equal(ref, as.double(out))
})

set.seed(200010) 
test_that("C++ rnorm works correctly", {
    vals <- scran:::test_rnorm(20000, 1, 1)
    expect_identical(length(vals), 20000L)
    expect_equal(mean(vals), 0, tol=0.01)
    expect_equal(var(vals), 1, tol=0.01)
    expect_identical(anyDuplicated(vals), 0L)

    vals2 <- scran:::test_rnorm(20000, 2, 1)
    expect_false(identical(vals, vals2))

    vals3 <- scran:::test_rnorm(20000, 1, 2)
    expect_false(identical(vals, vals3))
})

set.seed(20001)
test_that("correlateNull works with a design matrix", {
    # Constructing a reference function.
    REFFUN <- function(design, iters=1e3) {
        QR <- qr(design, LAPACK=TRUE)
        df <- nrow(design)-ncol(design)

        collected <- list()
        rand.state <- scran:::.setup_pcg_state(1e3)
        seeds <- rand.state$seeds[[1]]
        streams <- rand.state$streams[[1]]

        for (x in seq_along(seeds)) {
            vals <- scran:::test_rnorm(df*2L, seeds[[x]], streams[x])
            expect_identical(length(vals), df*2L)
    
            first.half <- qr.qy(QR, c(0,0, head(vals, df)))
            second.half <- qr.qy(QR, c(0, 0, tail(vals, df)))
            collected[[x]] <- cor(first.half, second.half, method="spearman")
        }
        sort(unlist(collected))
    }

    # A one-way layout.
    design <- model.matrix(~factor(rep(c(1,2), each=10)))
    
    set.seed(100)
    out1 <- REFFUN(design)
    set.seed(100)
    out2 <- correlateNull(design=design, iters=1e3)
    expect_equal(out1, as.double(out2))
    expect_equal(attr(out2, "design"), design)

    # A more complicated design.
    design <- model.matrix(~seq_len(10))
    set.seed(200)
    out1 <- REFFUN(design)
    set.seed(200)
    out2 <- correlateNull(design=design, iters=1e3)
    expect_equal(out1, as.double(out2))
    expect_equal(attr(out2, "design"), design)
})

test_that("correlateNull works with a blocking factor", {
    grouping <- rep(1:5, each=3)
    
    set.seed(100)
    out1 <- 0L
    for (gl in table(grouping)) { 
        out1 <- out1 + refnull(1e3, gl, resort=FALSE) * gl
    }
    out1 <- out1/length(grouping)
    out1 <- sort(out1)
    
    set.seed(100)
    out2 <- correlateNull(block=grouping, iters=1e3)
    expect_equal(out1, as.double(out2))
    expect_equal(attr(out2, "block"), grouping)
    expect_identical(attr(out2, "design"), NULL)
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
    expect_error(correlateNull(ncells=0), "number of cells should be greater than 2")
    expect_error(correlateNull(ncells=100, design=design), "cannot specify both 'ncells' and 'design'")
    expect_error(correlateNull(200, block=rep(1, 20)), "cannot specify")
    expect_error(correlateNull(200, design=cbind(rep(1, 20))), "cannot specify")
})
