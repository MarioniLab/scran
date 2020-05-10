# This tests the multiMarkerStats function.
# library(scran); library(testthat); source("test-multi-markers.R")

set.seed(1000)
library(scuttle)
sce <- mockSCE()
sce <- logNormCounts(sce)

# Any clustering method is okay, only using k-means for convenience.
kout <- kmeans(t(logcounts(sce)), centers=4)

tout <- findMarkers(sce, groups=kout$cluster, direction="up")
wout <- findMarkers(sce, groups=kout$cluster, direction="up", test="wilcox")
bout <- findMarkers(sce, groups=kout$cluster, direction="up", test="binom")

test_that("multiMarkerStats preserves single inputs correctly", {
    single <- multiMarkerStats(t=tout)

    for (i in seq_along(single)) {
        expect_equivalent(as.matrix(single[[i]][,1:3]), as.matrix(single[[i]][,1:3+3]))
        expect_equal(single[[i]][,1:3], tout[[i]][,1:3])

        lfc <- single[[i]][,-(1:6)]
        colnames(lfc) <- sub("^t.", "", colnames(lfc))
        expect_equal(lfc, tout[[i]][,-(1:3)])
    }
})

test_that("multiMarkerStats interleaves multiple inputs correctly", {
    combined <- multiMarkerStats(t=tout, wilcox=wout, binom=bout)

    for (i in seq_along(combined)) {
        curcom <- combined[[i]]
        expect_identical(
            colnames(curcom)[-(1:3)],
            as.character(rbind(
                paste0("t.", colnames(tout[[i]])),
                paste0("wilcox.", colnames(wout[[i]])),
                paste0("binom.", colnames(bout[[i]]))
            ))
        )

        curt <- tout[[i]][rownames(curcom),]
        curw <- wout[[i]][rownames(curcom),]
        curb <- bout[[i]][rownames(curcom),]

        expect_identical(curcom$Top, pmax(curt$Top, curw$Top, curb$Top))
        expect_identical(curcom$p.value, pmax(curt$p.value, curw$p.value, curb$p.value))
        
        expect_identical(curcom$t.Top, curt$Top)
        expect_identical(curcom$wilcox.Top, curw$Top)
        expect_identical(curcom$binom.Top, curb$Top)

        expect_equivalent(
            as.matrix(curt[,-(1:3)]), 
            as.matrix(curcom[,grep("^t\\..*logFC", colnames(curcom))])
        )
        expect_equivalent(
            as.matrix(curw[,-(1:3)]), 
            as.matrix(curcom[,grep("^wilcox\\..*AUC", colnames(curcom))])
        )
        expect_equivalent(
            as.matrix(curb[,-(1:3)]), 
            as.matrix(curcom[,grep("^binom\\..*logFC", colnames(curcom))])
        )
    }  
})

test_that("multiMarkerStats works correctly without 'Top' inputs", {
    tout2 <- findMarkers(sce, groups=kout$cluster, direction="up", pval.type="all")
    wout2 <- findMarkers(sce, groups=kout$cluster, direction="up", test="wilcox", pval.type="all")
    bout2 <- findMarkers(sce, groups=kout$cluster, direction="up", test="binom", pval.type="all")

    combined <- multiMarkerStats(t=tout2, wilcox=wout2, binom=bout2)
    for (i in seq_along(combined)) {
        curcom <- combined[[i]]
        expect_false("Top" %in% colnames(curcom))
        expect_false(is.unsorted(curcom$p.value))
    }
})

test_that("multiMarkerStats works correctly with log-transformed inputs", {
    ltout <- findMarkers(sce, groups=kout$cluster, log.p=TRUE, direction="up")
    lwout <- findMarkers(sce, groups=kout$cluster, log.p=TRUE, direction="up", test="wilcox")
    lbout <- findMarkers(sce, groups=kout$cluster, log.p=TRUE, direction="up", test="binom")

    ref <- multiMarkerStats(t=tout, wilcox=wout, binom=bout)
    combined <- multiMarkerStats(t=ltout, wilcox=lwout, binom=lbout)
    for (i in seq_along(combined)) {
        expect_equal(log(ref[[i]]$p.value), combined[[i]]$log.p.value)
        expect_equal(log(ref[[i]]$FDR), combined[[i]]$log.FDR)
    }
})

test_that("multiMarkerStats respects annotation correctly", { 
    stuff <- DataFrame(stuff=runif(nrow(sce)))
    rownames(stuff) <- rownames(sce)

    touta <- findMarkers(sce, groups=kout$cluster, direction="up", row.data=stuff)
    wouta <- findMarkers(sce, groups=kout$cluster, direction="up", test="wilcox", row.data=stuff)
    bouta <- findMarkers(sce, groups=kout$cluster, direction="up", test="binom", row.data=stuff)

    combined <- multiMarkerStats(t=touta, wilcox=wouta, binom=bouta, repeated="stuff")
    for (i in seq_along(combined)) {
        expect_equal(combined[[i]][rownames(sce),"stuff"], stuff$stuff)
    }
})

test_that("multiMarkerStats raises the expected set of errors", {
    toutx <- tout
    colnames(toutx[[1]]) <- paste0("A:", colnames(toutx[[1]]))
    expect_error(multiMarkerStats(t=toutx, wilcox=wout, binom=bout), "either all or no")

    toutx[[1]] <- toutx[[1]][,-1]
    expect_error(multiMarkerStats(t=toutx, wilcox=wout, binom=bout), "different numbers of columns")

    toutx[[1]] <- tout[[1]]
    rownames(toutx[[1]]) <- paste0("X", rownames(tout[[1]]))
    expect_error(multiMarkerStats(t=toutx, wilcox=wout, binom=bout), "row names")
    rownames(toutx[[1]]) <- NULL
    expect_error(multiMarkerStats(t=toutx, wilcox=wout, binom=bout), "non-NULL")
})
