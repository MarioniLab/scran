# Tests the pseudoBulkDGE function.
# library(testthat); library(scran); source("setup.R"); source("test-pseudo-dge.R")

set.seed(10000)
library(scater)
sce <- mockSCE(ncells=1000)
sce$samples <- gl(8, 125) # Pretending we have 8 samples.

# Making up some clusters.
sce <- logNormCounts(sce)
clusters <- kmeans(t(logcounts(sce)), centers=3)$cluster

# Creating a set of pseudo-bulk profiles:
info <- DataFrame(sample=sce$samples, cluster=clusters)
pseudo <- sumCountsAcrossCells(sce, info)

# Determining the experimental design for our 8 samples.
DRUG <- gl(2,4)
design <- model.matrix(~DRUG)
rownames(design) <- seq_len(8)

test_that("pseudoBulkDGE works correctly in vanilla cases", {
    # Spiking in DE for all clusters.
    pseudo2 <- pseudo
    xDRUG <- DRUG[as.integer(pseudo$sample)]
    assay(pseudo2)[1,xDRUG==1] <- assay(pseudo2)[1,xDRUG==1] * 100
    assay(pseudo2)[2,xDRUG==2] <- assay(pseudo2)[2,xDRUG==2] * 100

    out <- pseudoBulkDGE(pseudo2, 
       sample=pseudo2$sample, 
       label=pseudo2$cluster,
       design=design
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
    xDRUG <- DRUG[as.integer(pseudo$sample)]
    affected <- xDRUG==1 & pseudo2$cluster==2
    assay(pseudo2)[1,] <- 20 + as.integer(affected) * 100

    out <- pseudoBulkDGE(pseudo2,
       sample=pseudo2$sample,
       label=pseudo2$cluster,
       design=design
    )

    expect_true(out[[1]]$PValue[1] > 0.5)
    expect_true(out[[2]]$PValue[1] < 0.05)
    expect_true(out[[3]]$PValue[1] > 0.5)

    # Spits the dummy correctly.
    expect_error(
        pseudoBulkDGE(pseudo2,
            sample=pseudo2$sample,
            label=pseudo2$cluster,
            design=design[1:2,]
        ), "should have the same"
    )
})

test_that("pseudoBulkDGE handles the gene filtering correctly", {
    pseudo2 <- pseudo
    discard <- sample(nrow(pseudo), 50)
    assay(pseudo2)[discard,] <- 0

    out <- pseudoBulkDGE(pseudo2, 
        sample=pseudo2$sample, 
        label=pseudo2$cluster,
        design=design
    )

    ref <- pseudoBulkDGE(pseudo2[-discard,], 
        sample=pseudo2$sample, 
        label=pseudo2$cluster,
        design=design
    )

    for (i in names(out)) {
        expect_identical(out[[i]][-discard,], ref[[i]])
    }

    # Setting the group works.
    out2 <- pseudoBulkDGE(pseudo2,
        sample=pseudo$sample, 
        label=pseudo$cluster,
        design=design,
        condition=DRUG
    )

    expect_identical(out, out2)

    # Continues to work in the absence of row names.
    rownames(pseudo2) <- NULL

    out <- pseudoBulkDGE(pseudo2, 
        sample=pseudo2$sample, 
        label=pseudo2$cluster,
        design=design
    )

    for (i in names(out)) {
        curref <- ref[[i]]
        rownames(curref) <- NULL
        expect_identical(out[[i]][-discard,], curref)
    }
})

test_that("pseudoBulkDGE gracefully handles impossible comparisons", {
    discard <- pseudo$cluster == 3 & pseudo$sample %in% c(2:4, 6:8)
    pseudo2 <- pseudo[,!discard]

    out <- pseudoBulkDGE(pseudo2, 
        sample=pseudo2$sample, 
        label=pseudo2$cluster,
        design=design
    )

    expect_false(all(is.na(out[["3"]]$logFC)))
    expect_true(all(is.na(out[["3"]]$PValue)))

    # Checking that each cluster is truly independent of the others.
    ref <- pseudoBulkDGE(pseudo, 
        sample=pseudo$sample, 
        label=pseudo$cluster,
        design=design
    )
    
    expect_identical(ref[["1"]], out[["1"]])
    expect_identical(ref[["2"]], out[["2"]])

    # Smoothly handles situations where the comparison is impossible. 
    discard <- pseudo$cluster == 3 & pseudo$sample %in% c(5:8)
    pseudo2 <- pseudo[,!discard]

    out <- pseudoBulkDGE(pseudo2, 
        sample=pseudo2$sample, 
        label=pseudo2$cluster,
        design=design
    )

    expect_true(all(is.na(out[["3"]]$logFC)))
    expect_true(all(is.na(out[["3"]]$PValue)))
})

test_that("pseudoBulkDGE works with a log-fold change threshold", {
    out <- pseudoBulkDGE(pseudo, 
       sample=pseudo$sample, 
       label=pseudo$cluster,
       design=design,
       lfc=1
    )

    ref <- pseudoBulkDGE(pseudo, 
       sample=pseudo$sample, 
       label=pseudo$cluster,
       design=design
    )

    expect_false(identical(out, ref))

    # The lfc setting takes effect upon graceful failure.
    discard <- pseudo$cluster == 3 & pseudo$sample %in% c(2:4, 6:8)
    pseudo2 <- pseudo[,!discard]

    out <- pseudoBulkDGE(pseudo2,
        sample=pseudo2$sample,
        label=pseudo2$cluster,
        design=design, lfc=1,
    )

    expect_false(all(is.na(out[["3"]]$logFC)))
    expect_false(is.null(out[["3"]]$unshrunk.logFC))
    expect_false(all(is.na(out[["3"]]$unshrunk.logFC)))

    discard <- pseudo$cluster == 3 & pseudo$sample %in% c(5:8)
    pseudo2 <- pseudo[,!discard]

    out <- pseudoBulkDGE(pseudo2,
        sample=pseudo2$sample,
        label=pseudo2$cluster,
        design=design, lfc=1
    )

    expect_true(all(is.na(out[["3"]]$logFC)))
    expect_true(all(is.na(out[["3"]]$unshrunk.logFC)))
})

test_that("decideTestsPerLabel works correctly", {
    out <- pseudoBulkDGE(pseudo,
       sample=pseudo$sample, 
       label=pseudo$cluster,
       design=design
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
})
