# This tests the findMarkers function.
# require(scran); require(testthat); source("test-markers.R")

set.seed(70000000)
ncells <- 200
ngenes <- 250
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

library(scater)
rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
X <- logNormCounts(X)

test_that("findMarkers works correctly with subsetting and spikes", {   
    # Works with an SCE object.
    clust <- kmeans(t(exprs(X)), centers=3)
    out <- findMarkers(X, groups=clust$cluster)
    out2 <- findMarkers(exprs(X), groups=clust$cluster)
    expect_identical(out, out2)

    # Works with subsetting.
    out <- findMarkers(X, groups=clust$cluster, subset.row=100:1)
    out2 <- findMarkers(X[100:1,], groups=clust$cluster)
    expect_identical(out, out2)

    out <- findMarkers(X, groups=clust$cluster, subset.row=1:10, full.stats=TRUE) 
    out2 <- findMarkers(X[1:10,], groups=clust$cluster, full.stats=TRUE)
    expect_identical(out, out2)
    expect_identical(rownames(out[[1]]$stats.2), rownames(out[[1]])) # names propagate to internal DFs.

    # Switches to using 'subset.row' as the row names. 
    out <- findMarkers(X, groups=clust$cluster, subset.row=2:11, full.stats=TRUE, gene.names=NULL) 
    out2 <- findMarkers(unname(X)[2:11,], groups=clust$cluster, full.stats=TRUE, gene.names=2:11)
    expect_identical(out, out2)
    expect_identical(rownames(out[[1]]$stats.2), rownames(out[[1]])) 
    expect_identical(sort(as.integer(rownames(out[[1]]))), 2:11)

    # Repeating with a design matrix, to check that subsetting works in both branches for coefficient calculation.
    block <- factor(sample(2, ncol(X), replace=TRUE))
    design <- model.matrix(~block)[,-1,drop=FALSE]
    out.des <- findMarkers(exprs(X), groups=clust$cluster, design=design, subset.row=100:1)
    out.des2 <- findMarkers(exprs(X)[100:1,,drop=FALSE], groups=clust$cluster, design=design)
    expect_identical(out.des, out.des2)
})

test_that("findMarkers works correctly with row metadata", {
    clust <- kmeans(t(exprs(X)), centers=3)
    meta <- DataFrame(Y=runif(nrow(X)), row.names=rownames(dummy))
    out <- findMarkers(dummy, groups=clust$cluster, row.data=meta)

    for (i in seq_along(out)) {
        x <- out[[i]]
        expect_identical(x$Y, meta[rownames(x),"Y"])
    }

    # Handles it without names.
    out <- findMarkers(dummy, groups=clust$cluster, sorted=FALSE, gene.names=NULL, row.data=meta)
    for (i in seq_along(out)) {
        x <- out[[i]]
        expect_identical(rownames(x), as.character(seq_len(nrow(dummy))))
        expect_identical(x$Y, meta$Y)
    } 
})
