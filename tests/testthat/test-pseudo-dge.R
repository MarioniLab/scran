# Tests the pseudoBulkDGE function.
# library(testthat); library(scran); source("setup.R"); source("test-pseudo-dge.R")

set.seed(10000)
library(scuttle)
sce <- mockSCE(ncells=1000)
sce$samples <- gl(8, 125) # Pretending we have 8 samples.

# Making up some clusters.
sce <- logNormCounts(sce)
clusters <- kmeans(t(logcounts(sce)), centers=3)$cluster

# Creating a set of pseudo-bulk profiles:
info <- DataFrame(sample=sce$samples, cluster=clusters)
pseudo <- sumCountsAcrossCells(sce, info)
pseudo$DRUG <- gl(2,4)[pseudo$sample]

test_that("pseudoBulkDGE works correctly in vanilla cases", {
    # Spiking in DE for all clusters.
    pseudo2 <- pseudo
    xDRUG <- pseudo$DRUG
    assay(pseudo2)[1,xDRUG==1] <- assay(pseudo2)[1,xDRUG==1] * 100
    assay(pseudo2)[2,xDRUG==2] <- assay(pseudo2)[2,xDRUG==2] * 100

    out <- pseudoBulkDGE(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    expect_identical(names(out), as.character(sort(unique(clusters))))

    for (x in names(out)) {
       expect_identical(rownames(out[[x]]), rownames(sce))
       expect_true(out[[x]]$PValue[1] < 0.01)
       expect_true(out[[x]]$PValue[2] < 0.01)
       expect_true(out[[x]]$logFC[1] < -3)
       expect_true(out[[x]]$logFC[2] > 3)
    }

    # Spiking in DE for just one cluster.
    pseudo2 <- pseudo
    xDRUG <- pseudo$DRUG
    affected <- xDRUG==1 & pseudo2$cluster==2
    assay(pseudo2)[1,] <- 20 + as.integer(affected) * 100

    out <- pseudoBulkDGE(pseudo2,
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    expect_true(out[[1]]$PValue[1] > 0.5)
    expect_true(out[[2]]$PValue[1] < 0.05)
    expect_true(out[[3]]$PValue[1] > 0.5)

    # Handles design matrices as functions.
    out2 <- pseudoBulkDGE(pseudo2, label=pseudo2$cluster, 
        design=function(x) model.matrix(~DRUG, x), coef="DRUG2")

    expect_identical(out, out2)

    # Spits the dummy correctly.
    out <- pseudoBulkDGE(pseudo2, label=pseudo2$cluster, design=~BLAH)
    expect_identical(metadata(out)$failed, c("1", "2", "3"))
})

test_that("pseudoBulkDGE handles the gene filtering correctly", {
    pseudo2 <- pseudo
    discard <- sample(nrow(pseudo), 50)
    assay(pseudo2)[discard,] <- 0

    out <- pseudoBulkDGE(pseudo2, 
        label=pseudo2$cluster,
        design=~DRUG,
        coef="DRUG2"
    )

    ref <- pseudoBulkDGE(pseudo2[-discard,], 
        label=pseudo2$cluster,
        design=~DRUG,
        coef="DRUG2"
    )

    for (i in names(out)) {
        expect_identical(out[[i]][-discard,], ref[[i]])
    }

    # Setting the group works.
    out2 <- pseudoBulkDGE(pseudo2,
        label=pseudo2$cluster,
        design=~DRUG,
        coef="DRUG2",
        condition=pseudo$DRUG
    )

    expect_identical(out, out2)

    # Continues to work in the absence of row names.
    rownames(pseudo2) <- NULL

    out <- pseudoBulkDGE(pseudo2, 
        label=pseudo2$cluster,
        design=~DRUG,
        coef="DRUG2",
        include.intermediates=FALSE # because edgeR sticks row names onto everything.
    )

    for (i in names(out)) {
        curref <- ref[[i]]
        rownames(curref) <- NULL
        metadata(curref)$y <- NULL
        metadata(curref)$fit <- NULL
        expect_identical(out[[i]][-discard,], curref)
    }
})

test_that("pseudoBulkDGE gracefully handles impossible comparisons", {
    discard <- pseudo$cluster == 3 & pseudo$sample %in% c(2:4, 6:8)
    pseudo2 <- pseudo[,!discard]

    out <- pseudoBulkDGE(pseudo2, 
        label=pseudo2$cluster,
        design=~DRUG,
        coef="DRUG2"
    )

    expect_identical(metadata(out)$failed, "3")

    # Checking that each cluster's failure is truly independent of the others.
    ref <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2"
    )
    
    expect_identical(ref[["1"]], out[["1"]])
    expect_identical(ref[["2"]], out[["2"]])

    # Smoothly handles situations where the comparison is impossible. 
    discard <- pseudo$cluster == 3 & pseudo$sample %in% c(5:8)
    pseudo2 <- pseudo[,!discard]

    out <- pseudoBulkDGE(pseudo2, 
        label=pseudo2$cluster,
        design=~DRUG,
        coef="DRUG2"
    )

    expect_identical(metadata(out)$failed, "3")
})

test_that("pseudoBulkDGE works with a log-fold change threshold", {
    out <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
        lfc=1
    )

    ref <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2"
    )

    expect_false(identical(out, ref))
})

test_that("pseudoBulkDGE works with all the limma settings", {
    ref <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
        method="voom"
    )
    expect_identical(length(ref), 3L)

    out <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
        lfc=1,
        method="voom"
    )
    expect_identical(length(out), 3L)
    expect_false(identical(ref, out))

    out <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
        qualities=FALSE,
        method="voom"
    )
    expect_identical(length(out), 3L)
    expect_false(identical(ref, out))
})

test_that("contrast specification in the pseudoBulkDGE works", {
    ref <- pseudoBulkDGE(pseudo[1:10,], 
        label=pseudo$cluster,
        design=~DRUG,
        coef=2
    )

    out1 <- pseudoBulkDGE(pseudo[1:10,], 
        label=pseudo$cluster,
        design=~DRUG,
        contrast=c(0, 1)
    )

    out2 <- pseudoBulkDGE(pseudo[1:10,], 
        label=pseudo$cluster,
        design=~0 + DRUG,
        contrast=c(-1, 1)
    )

    out3 <- pseudoBulkDGE(pseudo[1:10,],
        label=pseudo$cluster,
        design=~0 + DRUG,
        contrast="DRUG2 - DRUG1"
    )

    for (i in seq_along(ref)) {
        expect_equal(ref[[i]]$logFC, out1[[i]]$logFC)
        expect_equal(ref[[i]]$logFC, out2[[i]]$logFC)
        expect_equal(ref[[i]]$logFC, out3[[i]]$logFC)
        expect_equal(ref[[i]]$PValue, out1[[i]]$PValue)
        expect_equal(ref[[i]]$PValue, out2[[i]]$PValue, tol=1e-6) # Need more generous tolerances on windows, who knows why.
        expect_equal(ref[[i]]$PValue, out3[[i]]$PValue, tol=1e-6)
    }

    # Voom is a bit different due to the approximation of the weights.
    ref <- pseudoBulkDGE(pseudo[1:10,], 
        label=pseudo$cluster,
        design=~0 + DRUG,
        method="voom",
        contrast=c(-1, 1)
    )

    out <- pseudoBulkDGE(pseudo[1:10,],
        label=pseudo$cluster,
        design=~0 + DRUG,
        method="voom",
        contrast="DRUG2 - DRUG1"
    )

    for (i in seq_along(ref)) {
        expect_equal(ref[[i]]$logFC, out[[i]]$logFC)
        expect_equal(ref[[i]]$PValue, out[[i]]$PValue)
    }
})

test_that("sorting in the pseudoBulkDGE works", {
    ref <- pseudoBulkDGE(pseudo, 
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
        sorted=TRUE,
        row.data=DataFrame(X=rownames(pseudo))
    )
   
    for (i in seq_len(3)) {
        expect_false(is.unsorted(ref[[i]]$PValue))
        expect_true("X" %in% colnames(ref[[i]]))
        expect_true(!is.null(metadata(ref[[i]])$y)) # metadata is preserved by the cbind.
    }
})

test_that("decideTestsPerLabel works correctly", {
    out <- pseudoBulkDGE(pseudo,
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
    )

    # Adding some DE genes to spice things up.
    for (i in names(out)) {
        chosen <- sample(nrow(out[[i]]), 50)
        out[[i]]$PValue[chosen] <- 0
    }

    dt <- decideTestsPerLabel(out)
    dt0 <- decideTestsPerLabel(out, lfc.field=NULL)
    expect_identical(abs(dt), dt0)

    dtg <- decideTestsPerLabel(out, method="global")
    expect_identical(dimnames(dtg), dimnames(dt))

    for (i in names(out)) {
        out[[i]]$logFC <- NULL
    }
    dt02 <- decideTestsPerLabel(out)
    expect_identical(dt02, dt0)

    # Works automatically with voom.
    out <- pseudoBulkDGE(pseudo,
        label=pseudo$cluster,
        design=~DRUG,
        coef="DRUG2",
        method="voom"
    )

    dt <- decideTestsPerLabel(out)
    dtp <- decideTestsPerLabel(out, pval.field="P.Value")
    expect_identical(dt, dtp)

    colnames(out[[1]]) <- head(LETTERS, ncol(out[[1]]))
    expect_error(decideTestsPerLabel(out), "pval.field")
})
