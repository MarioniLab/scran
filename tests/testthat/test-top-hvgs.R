# This tests the getTopHVGs function.
# library(testthat); library(scran); source("test-top-hvgs.R")

library(scuttle)
sce <- mockSCE()
sce <- logNormCounts(sce)

expect_identical_sorted <- function(x, y) expect_identical(sort(x), sort(y))

test_that("getTopHVGs works correctly", {
    stats <- modelGeneVar(sce)

    expect_identical_sorted(getTopHVGs(stats), 
        rownames(stats)[stats$bio > 0])

    expect_identical_sorted(getTopHVGs(stats, fdr.threshold=0.05), 
        rownames(stats)[stats$bio > 0 & stats$FDR <= 0.05])

    # Handles top choices correctly.
    expect_identical(getTopHVGs(stats, n=200, var.threshold=NULL), 
        head(rownames(stats)[order(-stats$bio)], 200))

    expect_identical(getTopHVGs(stats, n=200, prop=0.1, var.threshold=NULL), 
        head(rownames(stats)[order(-stats$bio)], 0.1*nrow(stats)))

    expect_identical(getTopHVGs(stats, n=2000, prop=0.001, var.threshold=NULL), 
        head(rownames(stats)[order(-stats$bio)], 2000))

    expect_identical(getTopHVGs(stats, n=NULL, prop=0.1, var.threshold=NULL), 
        head(rownames(stats)[order(-stats$bio)], 0.1*nrow(stats)))

    expect_identical(getTopHVGs(stats, n=Inf, var.threshold=NULL), 
        rownames(stats)[order(-stats$bio)])
})

test_that("getTopHVGs handles unnamed inputs correctly", {
    stats <- modelGeneVar(sce)
    rownames(stats) <- NULL

    expect_identical_sorted(getTopHVGs(stats), which(stats$bio > 0))

    expect_identical_sorted(getTopHVGs(stats, fdr.threshold=0.05),
        which(stats$bio > 0 & stats$FDR <= 0.05))

    expect_identical(getTopHVGs(stats, var.threshold=NULL),
        head(order(stats$bio, decreasing=TRUE), 2000))
})

test_that("getTopHVGs handles NA values", {
    stats <- modelGeneVar(sce)
    stats$bio[1] <- 10000
    stats$FDR[1] <- NA
    stats$bio[2] <- NA
    stats$FDR[2] <- 0

    expect_identical_sorted(
        getTopHVGs(stats, fdr.threshold=0.05, var.threshold=NULL),
        getTopHVGs(stats[-1,], fdr.threshold=0.05, var.threshold=NULL)
    )

    expect_identical_sorted(
        getTopHVGs(stats),
        getTopHVGs(stats[-2,])
    )

    expect_identical_sorted(
        getTopHVGs(stats, fdr.threshold=0.05),
        getTopHVGs(stats[-(1:2),], fdr.threshold=0.05)
    )
})

test_that("getTopHVGs works correctly with CV2", {
    stats2 <- modelGeneCV2(sce)
    
    expect_identical_sorted(getTopHVGs(stats2, var.field="ratio", var.threshold=1),
        rownames(stats2)[stats2$ratio > 1])
    
    expect_identical_sorted(getTopHVGs(stats2, var.field="ratio", var.threshold=1, fdr.threshold=0.05),
        rownames(stats2)[stats2$ratio > 1 & stats2$FDR <= 0.05])

    expect_identical_sorted(getTopHVGs(stats2, var.field="ratio", n=200, var.threshold=NULL), 
        head(rownames(stats2)[order(-stats2$ratio)], 200))

    expect_identical_sorted(getTopHVGs(stats2, var.field="ratio", n=Inf, var.threshold=NULL), 
        rownames(stats2)[order(-stats2$ratio)])
})
