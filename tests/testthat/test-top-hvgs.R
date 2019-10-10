# This tests the getTopHVGs function.
# library(testthat); library(scran); source("test-top-hvgs.R")

library(scater)
sce <- mockSCE()
sce <- logNormCounts(sce)

expect_identical_sorted <- function(x, y) expect_identical(sort(x), sort(y))

test_that("getTopHVGs works correctly", {
    stats <- modelGeneVar(sce)

    expect_identical_sorted(getTopHVGs(stats), 
        rownames(stats)[stats$bio > 0])

    expect_identical_sorted(getTopHVGs(stats, fdr.threshold=0.05), 
        rownames(stats)[stats$bio > 0 & stats$FDR <= 0.05])

    expect_identical_sorted(getTopHVGs(stats, n=200, var.threshold=NULL), 
        head(rownames(stats)[order(-stats$bio)], 200))

    expect_identical_sorted(getTopHVGs(stats, n=Inf, var.threshold=NULL), 
        rownames(stats)[order(-stats$bio)])
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
