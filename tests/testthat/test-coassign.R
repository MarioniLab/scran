# Tests the coassignProb() function.
# library(testthat); library(scran); source('setup.R'); source('test-coassign.R')

library(scater)
sce <- mockSCE(ncells=200)
sce <- logNormCounts(sce)

set.seed(8000)
clust1 <- kmeans(t(logcounts(sce)),3)$cluster
clust2 <- kmeans(t(logcounts(sce)),5)$cluster

test_that("coassignProb gives the expected output", {
    prob <- coassignProb(clust1, clust2)

    expect_identical(dim(prob), c(3L, 3L))
    expect_identical(rownames(prob), as.character(1:3))
    expect_identical(colnames(prob), as.character(1:3))

    expect_true(all(is.na(prob[lower.tri(prob)])))

    vals <- prob[upper.tri(prob, diag=TRUE)]
    expect_true(all(!is.na(vals)))
    expect_true(all(vals >= 0))
    expect_true(all(vals <= 1))

    # Correct values emitted.
    prob <- coassignProb(clust1, clust1)
    expect_true(all(diag(prob)==1))
    expect_true(all(prob[upper.tri(prob)]==0))
})

test_that("coassignProb handles silly inputs", {
    ref <- coassignProb(clust1, clust2)

    # Respects empty levels in 'ref'.
    clustX <- factor(clust1, levels=1:5)
    out <- coassignProb(clustX, clust2)

    expect_identical(dim(out), c(5L, 5L))
    expect_identical(rownames(out), as.character(1:5))
    expect_identical(colnames(out), as.character(1:5))
    expect_identical(out[1:3,1:3], ref)

    # Ignores empty levels in 'alt'.
    clustY <- factor(clust2, levels=1:5)
    out <- coassignProb(clust1, clustY)
    expect_identical(ref, out)

    # Behaves sensibly with zero-length inputs.
    out <- coassignProb(clust1[0], clust2[0])
    expect_identical(dim(out), c(0L, 0L))
})

test_that("coassignProb works with summarization", {
    ref <- coassignProb(clust1, clust2)
    stats <- coassignProb(clust1, clust2, summarize=TRUE)
    expect_equivalent(diag(ref), stats$self)

    ref0 <- ref
    diag(ref0) <- 0
    expect_equivalent(pmax(colMaxs(ref0, na.rm=TRUE), rowMaxs(ref0, na.rm=TRUE)), stats$other)

    # Works correctly when there are missing levels.
    clustX <- factor(clust1, levels=1:5)
    stats2 <- coassignProb(clustX, clust2, summarize=TRUE)
    expect_identical(nrow(stats2), 5L)
    expect_identical(stats2[1:3,], stats)
})
