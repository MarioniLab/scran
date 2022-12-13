# Tests the quickSubCluster utility.
# require(scran); require(testthat); source("setup.R"); source("test-subclust.R")

library(scran)

set.seed(30000)
ncells <- 700
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)

set.seed(30001)
test_that("quickSubCluster's defaults are consistent with quickCluster's defaults", {
    output <- quickCluster(dummy, min.size=0)
    output2 <- quickSubCluster(dummy, groups=rep(1, ncol(dummy)))
    expect_identical(as.character(output), sub(".*\\.", "", output2[[1]]$subcluster))

    sampling <- sample(3, ncells, replace=TRUE)
    output <- quickSubCluster(dummy, groups=sampling)
    for (i in unique(sampling)) {
        ref <- quickCluster(dummy[,i==sampling], min.size=0)
        expect_identical(as.character(ref), sub(".*\\.", "", output[[i]]$subcluster))
    }
})

set.seed(30001)
test_that("quickSubCluster's metadata output makes sense", {
    sampling <- sample(3, ncells, replace=TRUE)
    output <- quickSubCluster(dummy, groups=sampling)

    for (i in unique(sampling)) {
        expect_identical(metadata(output)$index[[i]], which(sampling==i))
        expect_identical(metadata(output)$subcluster[metadata(output)$index[[i]]], output[[i]]$subcluster)
    }

    expect_identical(sub("\\..*", "", metadata(output)$subcluster), as.character(sampling))

    raw.out <- quickSubCluster(dummy, groups=sampling, simplify=TRUE)
    expect_identical(metadata(output)$subcluster, raw.out)
})

set.seed(30002)
test_that("quickSubCluster behaves correctly upon changing the assay", {
    dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
    sampling <- sample(3, ncells, replace=TRUE)
    sce <- SingleCellExperiment(list(whee=dummy, blah=scuttle::normalizeCounts(dummy)))

    output <- quickSubCluster(dummy, groups=sampling)
    output2 <- quickSubCluster(sce, groups=sampling, assay.type="whee")
    expect_identical(lapply(output, reducedDim), lapply(output2, reducedDim))

    # Checking that we get the same output when we turn off normalization.
    expect_error(output <- quickSubCluster(scuttle::normalizeCounts(dummy), groups=sampling, normalize=FALSE), NA)
    output2 <- quickSubCluster(sce, groups=sampling, normalize=FALSE, assay.type="blah")
    expect_identical(lapply(output, reducedDim), lapply(output2, reducedDim))
})

set.seed(30003)
test_that("quickSubCluster avoids subclustering with too few cells", {
    ncells <- 99
    dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
    sampling <- sample(2, ncells, replace=TRUE)

    suppressWarnings(output <- quickSubCluster(dummy, groups=sampling))
    has.sub1 <- any(grepl("\\.", output[[1]]$subcluster))
    has.sub2 <- any(grepl("\\.", output[[2]]$subcluster))
    expect_true(has.sub1!=has.sub2)

    # Avoids crashing with one-cell clusters.
    test <- quickSubCluster(dummy, groups=rep(1:2, c(1, ncells-1)))
    expect_equivalent(counts(test[[1]]), dummy[,1,drop=FALSE])
})

set.seed(30001)
test_that("quickSubCluster restrictions work as expected", {
    sampling <- sample(LETTERS[1:3], ncells, replace=TRUE)

    set.seed(100)
    suppressWarnings(full <- quickSubCluster(dummy, groups=sampling))
    expect_identical(names(full), LETTERS[1:3])

    set.seed(100)
    suppressWarnings(res <- quickSubCluster(dummy, groups=sampling, restricted="A"))
    expect_identical(names(res), "A")
    expect_identical(res$A, full$A)

    idx <- which(sampling == "A")
    expect_identical(metadata(res)$index, list(A=idx))
    expect_identical(metadata(res)$subcluster[idx], metadata(full)$subcluster[idx])
    expect_identical(metadata(res)$subcluster[-idx], sampling[-idx])
})

