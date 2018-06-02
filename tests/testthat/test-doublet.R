# This tests the doublet-discovery machinery, 
# though there is little to do other than run the function.
# library(scran); library(testthat); source("test-doublet.R")

set.seed(9900001)
ngenes <- 100
mu1 <- 2^rexp(ngenes)
mu2 <- 2^rnorm(ngenes)

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes)
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes)
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes)

counts <- cbind(counts.1, counts.2, counts.m)
clusters <- rep(1:3, c(ncol(counts.1), ncol(counts.2), ncol(counts.m)))

test_that("doubletCluster works correctly with vanilla tests", {
	dbl <- doubletCluster(counts, clusters)
    expect_identical(rownames(dbl)[1], "3")
    expect_identical(dbl$source1[1], "2")
    expect_identical(dbl$source2[1], "1")
	
	# Checking the relative library sizes.
    ls1 <- median(colSums(counts.1))
    ls2 <- median(colSums(counts.2))
    ls3 <- median(colSums(counts.m))

	expect_equal(dbl$lib.size1[1], ls2/ls3)
	expect_equal(dbl$lib.size2[1], ls1/ls3)

    # Checking the proportions.
    expect_equal(dbl$prop,
        as.integer(table(clusters)[rownames(dbl)])/length(clusters))

    # Checking that p-values are reverse-sorted.
    expect_false(is.unsorted(-dbl$p.value))

    # Checking that we get equivalent results with character cluster input.
    re.clusters <- LETTERS[clusters]
    re.dbl <- doubletCluster(counts, re.clusters)
    dbl2 <- dbl
    rownames(dbl2) <- LETTERS[as.integer(rownames(dbl2))]
    dbl2$source1 <- LETTERS[as.integer(dbl2$source1)]
    dbl2$source2 <- LETTERS[as.integer(dbl2$source2)]
    expect_identical(dbl2, re.dbl)
})

test_that("doubletCluster works correctly with row subsets", {
    chosen <- sample(ngenes, 20)
    dbl0 <- doubletCluster(counts, clusters, subset.row=chosen)
    ref <- doubletCluster(counts[chosen,], clusters)
    ref$best <- as.character(chosen)[as.integer(ref$best)]
    expect_identical(dbl0, ref)

    # Trying out empty rows.
    out <- doubletCluster(counts[0,], clusters)
    expect_identical(nrow(out), nrow(ref))
    expect_true(all(is.na(out$best)))
    expect_true(all(is.na(out$p.value)))
    expect_true(all(out$N==0L))

    # While we're here, trying out empty columns.
    out <- doubletCluster(counts[,0], clusters[0])
    expect_identical(out, ref[0,])
})

test_that("doubletCluster works correctly with SingleCellExperiment", {
    sce <- SingleCellExperiment(list(counts=counts))
	ref <- doubletCluster(counts, clusters)
	dbl <- doubletCluster(sce, clusters)
    expect_identical(ref, dbl)

    # With a different assay.
    assay(sce, "whee") <- counts + rpois(length(counts), lambda=2)
	ref2 <- doubletCluster(assay(sce, "whee"), clusters)
	dbl2 <- doubletCluster(sce, clusters, assay.type="whee")
    expect_identical(ref2, dbl2)

    # With spike-ins that get used.
    isSpike(sce, "ERCC") <- sample(nrow(sce), 20)
	dbl3 <- doubletCluster(sce, clusters, get.spikes=TRUE)
    expect_identical(ref, dbl3)

    # ... or ignored.
	dbl4 <- doubletCluster(sce, clusters)
    ref3 <- doubletCluster(counts(sce), clusters, subset.row=!isSpike(sce))
    expect_identical(ref3, dbl4)

    # With both spike-ins _and_ subset.row specified.
    keep <- c(sample(which(isSpike(sce)), 10), sample(which(!isSpike(sce)), 10))
   	dbl5 <- doubletCluster(sce, clusters, subset.row=keep)
    ref4 <- doubletCluster(counts(sce), clusters, subset.row=setdiff(keep, which(isSpike(sce))))
    expect_identical(ref4, dbl5)
})

