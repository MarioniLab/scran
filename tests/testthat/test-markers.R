# This tests the findMarkers function.
# require(scran); require(testthat); source("test-markers.R")

set.seed(70000000)
ncells <- 200
ngenes <- 250
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

library(scuttle)
rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
X <- logNormCounts(X)

test_that("findMarkers dispatches correctly", {   
    clust <- kmeans(t(logcounts(X)), centers=3)
    out <- findMarkers(X, groups=clust$cluster)
    out2 <- findMarkers(logcounts(X), groups=clust$cluster)
    expect_identical(out, out2)

    Xbase <- as(X, "SummarizedExperiment")
    rownames(Xbase) <- rownames(X) # TODO: Bioconductor/SummarizedExperiment#29
    out3 <- findMarkers(Xbase, groups=clust$cluster)
    expect_identical(out, out3)

    Xlab <- X
    colLabels(Xlab) <- clust$cluster
    out4 <- findMarkers(Xlab)
    expect_identical(out, out4)
})

test_that("findMarkers works correctly with subsetting and spikes", {   
    clust <- kmeans(t(logcounts(X)), centers=3)

    # Works with subsetting.
    out <- findMarkers(X, groups=clust$cluster, subset.row=100:1)
    out2 <- findMarkers(X[100:1,], groups=clust$cluster)
    expect_identical(out, out2)

    out <- findMarkers(X, groups=clust$cluster, subset.row=1:10, full.stats=TRUE) 
    out2 <- findMarkers(X[1:10,], groups=clust$cluster, full.stats=TRUE)
    expect_identical(out, out2)
    expect_identical(rownames(out[[1]]$stats.2), rownames(out[[1]])) # names propagate to internal DFs.

    # Repeating with a design matrix, to check that subsetting works in both branches for coefficient calculation.
    block <- factor(sample(2, ncol(X), replace=TRUE))
    design <- model.matrix(~block)[,-1,drop=FALSE]
    out.des <- findMarkers(logcounts(X), groups=clust$cluster, design=design, subset.row=100:1)
    out.des2 <- findMarkers(logcounts(X)[100:1,,drop=FALSE], groups=clust$cluster, design=design)
    expect_identical(out.des, out.des2)
})

test_that("findMarkers works correctly with row metadata", {
    clust <- kmeans(t(logcounts(X)), centers=3)
    meta <- DataFrame(Y=runif(nrow(X)), row.names=rownames(dummy))
    out <- findMarkers(dummy, groups=clust$cluster, row.data=meta)

    for (i in seq_along(out)) {
        x <- out[[i]]
        expect_identical(x$Y, meta[rownames(x),"Y"])
    }

    # Handles it without sorting.
    out <- findMarkers(dummy, groups=clust$cluster, sorted=FALSE, row.data=meta)
    for (i in seq_along(out)) {
        x <- out[[i]]
        expect_identical(rownames(x), paste0("X", seq_len(nrow(dummy))))
        expect_identical(x$Y, meta$Y)
    } 

    # Errors out properly.
    nonames <- meta
    rownames(nonames) <- NULL
    expect_error(findMarkers(dummy, groups=clust$cluster, row.data=nonames), "inconsistent or NULL")
    expect_error(findMarkers(dummy, groups=clust$cluster, sorted=FALSE, row.data=nonames), "inconsistent")
})

test_that("findMarkers works correctly with row metadata as a list", {
    clust <- kmeans(t(logcounts(X)), centers=3)
    meta <- summaryMarkerStats(dummy, groups=clust$cluster)
    out <- findMarkers(dummy, groups=clust$cluster, sorted=FALSE, row.data=meta)

    for (i in seq_along(out)) {
        x <- out[[i]]
        curclust <- names(out)[i]
        expect_equal(x$self.average, rowMeans(dummy[,curclust==clust$cluster]))
        expect_equal(x$self.detected, rowMeans(dummy[,curclust==clust$cluster] != 0))
    }

    # Handles sorting properly.
    out2 <- findMarkers(dummy, groups=clust$cluster, sorted=TRUE, row.data=meta)
    for (i in seq_along(out)) {
        cur1 <- out[[i]]
        cur2 <- out2[[i]]
        expect_identical(cur1[rownames(cur2),], cur2)
    }

    # Auto-computes summaries.
    auto <- findMarkers(dummy, groups=clust$cluster, sorted=FALSE, add.summary=TRUE)
    expect_identical(auto, out)
})

test_that("findMarkers and getTopMarkers work correctly", {
    clust <- kmeans(t(logcounts(X)), centers=3)
    stats <- pairwiseTTests(dummy, groups=clust$cluster)

    out <- findMarkers(dummy, groups=clust$cluster)
    top <- getTopMarkers(stats[[1]], stats[[2]], pairwise=FALSE, fdr.threshold=NULL)
    ref <- lapply(out, FUN=function(x) rownames(x)[x$Top <= 10]) 
    expect_identical(as.list(top), ref)

    out <- findMarkers(dummy, groups=clust$cluster, pval.type="all")
    top <- getTopMarkers(stats[[1]], stats[[2]], pairwise=FALSE, pval.type="all", fdr.threshold=NULL)
    ref <- lapply(out, FUN=function(x) rownames(x)[1:10])
    expect_identical(as.list(top), ref)

    top <- getTopMarkers(stats[[1]], stats[[2]], pairwise=FALSE, pval.type="all")
    ref <- lapply(out, FUN=function(x) head(rownames(x)[x$FDR <= 0.05], 10))
    expect_identical(as.list(top), ref)

    # Checking with pairwise=TRUE.
    out <- getTopMarkers(stats[[1]], stats[[2]], pairwise=TRUE, fdr.threshold=NULL)
    expect_identical(unique(lengths(out)), 3L)
    expect_equivalent(do.call(cbind, lapply(out, lengths)), (1 - diag(3)) * 10L)
    expect_identical(unique(lapply(out, names)), list(as.character(1:3)))

    alt <- getTopMarkers(stats[[1]], stats[[2]], pairwise=FALSE, fdr.threshold=NULL)
    expect_identical(lapply(lapply(lapply(out, unlist), unique), sort), lapply(alt, sort))

    # Checking some genes get thrown out by the FDR filter.
    bounded <- getTopMarkers(stats[[1]], stats[[2]], pairwise=TRUE)
    expect_true(all(unlist(lapply(bounded, lengths)) <= unlist(lapply(out, lengths))))
    expect_true(length(unlist(bounded)) <= length(unlist(out)))

    # Checking it tolerates NAs.
    nastats <- stats
    nastats[[1]]$FDR[1] <- NA
    naive <- getTopMarkers(nastats[[1]], nastats[[2]], pairwise=TRUE)

    nastats[[1]]$FDR[1] <- 1
    ref <- getTopMarkers(nastats[[1]], nastats[[2]], pairwise=TRUE)
    expect_identical(naive, ref)
})

test_that("findMarkers and getMarkerEffects work correctly", {
    clust <- kmeans(t(logcounts(X)), centers=3)

    out <- findMarkers(dummy, groups=clust$cluster)
    eff <- getMarkerEffects(out[[1]])
    expect_type(eff, "double")
    expect_identical(colnames(eff), as.character(2:3))

    # Removes NAs properly.
    copy <- out[[1]]
    ref <- getMarkerEffects(copy)
    copy$logFC.2 <- NA
    eff <- getMarkerEffects(copy, remove.na.col=TRUE)
    expect_identical(ref[,-1,drop=FALSE], eff)

    # Works for Wilcox tests.
    out <- findMarkers(dummy, groups=clust$cluster, test.type="wilcox")
    eff <- getMarkerEffects(out[[2]], prefix="AUC")
    expect_type(eff, "double")
    expect_identical(colnames(eff), as.character(c(1,3)))
})
