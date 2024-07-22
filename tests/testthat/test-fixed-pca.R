# Tests the fixedPCA function.
# library(testthat); library(scran); source("test-fixed-pca.R")

library(scuttle)
sce <- mockSCE()
sce <- logNormCounts(sce)
library(BiocSingular)

test_that("fixedPCA works correctly", {
    set.seed(100)
    sce2 <- fixedPCA(sce, subset.row=NULL)
    set.seed(100)
    ref <- runPCA(t(logcounts(sce2)), rank=50, BSPARAM=bsparam())
    expect_equal(unclass(reducedDim(sce2))[,], ref$x)

    set.seed(100)
    sce2 <- fixedPCA(sce, subset.row=1:200)
    set.seed(100)
    ref <- runPCA(t(logcounts(sce2)[1:200,]), rank=50, BSPARAM=bsparam())
    expect_equal(unclass(reducedDim(sce2))[,], ref$x)
    expect_equal(logcounts(sce), logcounts(sce2))

    # Doesn't preserve shape if we don't ask.
    set.seed(100)
    sce2 <- fixedPCA(sce, subset.row=1:200, preserve.shape=FALSE)
    expect_equal(unclass(reducedDim(sce2))[,], ref$x)
    expect_equal(logcounts(sce2), logcounts(sce)[1:200,])

    set.seed(100)
    sce <- fixedPCA(sce, rank=20, subset.row=1:50)
    set.seed(100)
    ref <- runPCA(t(logcounts(sce)[1:50,]), rank=20, BSPARAM=bsparam())
    expect_equal(unclass(reducedDim(sce))[,], ref$x)
})

test_that("fixedPCA works correctly with low rank approximations", {
    set.seed(100)
    sce <- fixedPCA(sce, subset.row=NULL)
    set.seed(100)
    sce2 <- fixedPCA(sce, value="lowrank", subset.row=NULL)
    rot <- attr(reducedDim(sce), "rotation")
    expect_identical(as.matrix(assay(sce2, "lowrank")[1:10,]), tcrossprod(rot[1:10,], reducedDim(sce)))

    # Works with subsetting.
    set.seed(100)
    sce3 <- fixedPCA(rbind(sce, sce[1:10,]), subset.row=seq_len(nrow(sce)), value="lowrank")
    expect_identical(assay(sce2, "lowrank"), assay(sce3, "lowrank")[seq_len(nrow(sce)),])
    expect_equal(assay(sce2, "lowrank")[1:10,], assay(sce3, "lowrank")[nrow(sce)+1:10,], tol=1e-6)

    #  Won't preserve the shape.
    set.seed(100)
    sce4 <- fixedPCA(rbind(sce, sce[1:10,]), subset.row=seq_len(nrow(sce)), value="lowrank", preserve.shape=FALSE)
    expect_identical(assay(sce2, "lowrank"), assay(sce4, "lowrank"))
})

test_that("fixedPCA warns when subset.row is not specified", {
    expect_warning(fixedPCA(sce), "subset.row")
})

