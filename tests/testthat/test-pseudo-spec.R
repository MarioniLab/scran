# This tests the pseudoBulkSpecific functions.
# library(testthat); library(scran); source("test-pseudo-spec.R")

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

test_that("pseudoBulkSpecific works correctly in vanilla cases", {
    # Spiking in DE for all clusters.
    pseudo2 <- pseudo
    xDRUG <- pseudo$DRUG
    assay(pseudo2)[1,xDRUG==1] <- assay(pseudo2)[1,xDRUG==1] * 100
    assay(pseudo2)[2,xDRUG==2] <- assay(pseudo2)[2,xDRUG==2] * 100

    ref <- pseudoBulkDGE(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    out <- pseudoBulkSpecific(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    expect_identical(names(ref), names(out))
    for (i in names(ref)) {
        left <- ref[[i]]
        right <- out[[i]]
        expect_true(all(left$PValue <= right$PValue))
        expect_identical(left$LogFC, right$LogFC)
        expect_true(all(right$PValue[1:2] > 0.01))
        expect_true(all(left$PValue[1:2] < 0.01))
    }

    # Also works for voom.
    ref <- pseudoBulkDGE(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2",
       method="voom"
    )

    out <- pseudoBulkSpecific(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2",
       method="voom"
    )

    expect_identical(names(ref), names(out))
    for (i in names(ref)) {
        left <- ref[[i]]
        right <- out[[i]]
        expect_true(all(left$P.Value <= right$P.Value))
        expect_identical(left$LogFC, right$LogFC)
        expect_true(all(right$P.Value[1:2] > 0.05))
        expect_true(all(left$P.Value[1:2] < 0.05))
    }
})

test_that("pseudoBulkSpecific gives the same answers with a reference", {
    out <- pseudoBulkSpecific(pseudo, 
       label=pseudo$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    ref <- pseudoBulkDGE(pseudo, 
       label=pseudo$cluster,
       design=~DRUG,
       coef="DRUG2"
    )
    metadata(ref)$tag <- "I'm here!"
    metadata(ref[[1]])$tag <- "I'm still here!"

    out2 <- pseudoBulkSpecific(pseudo,
       label=pseudo$cluster,
       design=~DRUG,
       coef="DRUG2",
       reference=ref
    )

    expect_identical(metadata(out2)$tag, "I'm here!")
    expect_identical(metadata(out2[[1]])$tag, "I'm still here!")

    metadata(out2)$tag <- NULL
    metadata(out2[[1]])$tag <- NULL
    expect_identical(out, out2)
})

test_that("pseudoBulkSpecific works correctly with sorting", {
    ref <- pseudoBulkDGE(pseudo, 
       label=pseudo$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    out <- pseudoBulkSpecific(pseudo, 
       label=pseudo$cluster,
       design=~DRUG,
       coef="DRUG2",
       reference=ref
    )

    out2 <- pseudoBulkSpecific(pseudo,
       label=pseudo$cluster,
       design=~DRUG,
       coef="DRUG2",
       reference=ref,
       sorted=TRUE
    )

    for (i in names(ref)) {
        left <- out[[i]]
        right <- out2[[i]]
        expect_identical(left[order(left$PValue),], right)
    }
})

test_that("pseudoBulkSpecific works correctly with zero replacement", {
    pseudo2 <- pseudo
    assay(pseudo2)[1,pseudo2$cluster!=1] <- 0

    ref <- pseudoBulkDGE(pseudo2,
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2"
    )

    out <- pseudoBulkSpecific(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2",
       reference=ref,
    )

    out2 <- pseudoBulkSpecific(pseudo2, 
       label=pseudo2$cluster,
       design=~DRUG,
       coef="DRUG2",
       reference=ref,
       missing.as.zero=TRUE
    )

    expect_identical(out[[1]]$OtherAverage[1], NA_real_)
    expect_identical(out2[[1]]$OtherAverage[1], 0)
})

test_that("pseudoBulkSpecific errors out correctly", {
    expect_error(pseudoBulkSpecific(pseudo, 
       label=pseudo$cluster,
       design=~DRUG,
       coef=1:2
    ), "cannot be specified")
})
