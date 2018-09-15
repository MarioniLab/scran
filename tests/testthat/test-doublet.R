# This tests the doublet discovery machinery in scran.
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

RENAMER <- function(val, fields, mapping) 
# A convenience function for remapping internal fields upon subsetting or renaming.
# This is necessary for some equality checks below.
{
    new.pairs <- val$all.pairs
    for (f in fields) {
        val[[f]] <- mapping[as.integer(val[[f]])]
        for (i in seq_along(new.pairs)) {
            new.pairs[[i]][[f]] <- mapping[as.integer(new.pairs[[i]][[f]])]
        }
    }
    val$all.pairs <- new.pairs
    val
}

###################################################

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
    expect_equal(dbl$prop, as.integer(table(clusters)[rownames(dbl)])/length(clusters))

    # Checking that p-values are reverse-sorted.
    expect_false(is.unsorted(-dbl$p.value))

    # Checking that we get equivalent results with character cluster input.
    re.clusters <- LETTERS[clusters]
    re.dbl <- doubletCluster(counts, re.clusters)

    dbl2  <- RENAMER(dbl, c("source1", "source2"), LETTERS)
    rownames(dbl2) <- LETTERS[as.integer(rownames(dbl2))]
    names(dbl2$all.pairs) <- LETTERS[as.integer(names(dbl2$all.pairs))]
    expect_identical(dbl2, re.dbl)
})

test_that("doubletCluster agrees with a reference implementation", {
    mu3 <- 2^rnorm(ngenes)
    counts.3 <- matrix(rpois(ngenes*100, mu3), nrow=ngenes)
    counts <- cbind(counts.1, counts.2, counts.3, counts.m)
    clusters <- rep(1:4, c(ncol(counts.1), ncol(counts.2), ncol(counts.3), ncol(counts.m)))

    dbl <- doubletCluster(counts, clusters)
    ref <- findMarkers(scater::normalizeCounts(counts, scater::librarySizeFactors(counts)), clusters, full.stats=TRUE)

    for (x in rownames(dbl)) {
        stats <- ref[[x]]
        all.pops <- setdiff(rownames(dbl), x)
        combos <- combn(all.pops, 2)

        # Effectively a re-implentation of the two inner loops.
        collected <- apply(combos, 2, function(chosen) {
            fields <- paste0("stats.", chosen)
            stats1 <- stats[[fields[1]]]
            stats2 <- stats[[fields[2]]]
            p <- pmax(exp(stats1$log.p.value), exp(stats2$log.p.value))
            p[sign(stats1$logFC)!=sign(stats2$logFC)] <- 1
            adj.p <- p.adjust(p, method="BH")
            data.frame(best=rownames(stats)[which.min(p)], p.val=min(adj.p), N=sum(adj.p <= 0.05), stringsAsFactors=FALSE)
        })

        collected <- do.call(rbind, collected)
        o <- order(collected$N)

        obs <- dbl[x,"all.pairs"][[1]]
        expect_identical(obs$source1, pmax(combos[2,], combos[1,])[o])
        expect_identical(obs$source2, pmin(combos[1,], combos[2,])[o])
        expect_identical(obs$N, collected$N[o])
        expect_identical(obs$best, collected$best[o])
        expect_equal(obs$p.value, collected$p.val[o])

        to.use <- o[1]
        expect_identical(dbl[x,"N"], collected[to.use, "N"])
        expect_equal(dbl[x,"p.value"], collected[to.use, "p.val"])
        expect_identical(dbl[x,"best"], collected[to.use, "best"])
        expect_identical(sort(c(dbl[x,"source1"],dbl[x,"source2"])), sort(combos[,to.use]))
    }
})

test_that("doubletCluster works correctly with row subsets", {
    chosen <- sample(ngenes, 20)
    dbl0 <- doubletCluster(counts, clusters, subset.row=chosen)
    ref <- doubletCluster(counts[chosen,], clusters)
    ref <- RENAMER(ref, "best", as.character(chosen))
    expect_identical(dbl0, ref)

    # Trying out empty rows.
    out <- doubletCluster(counts[0,], clusters)
    expect_identical(nrow(out), nrow(ref))
    expect_true(all(is.na(out$best)))
    expect_true(all(is.na(out$p.value)))
    expect_true(all(out$N==0L))

    # While we're here, trying out empty columns.
    expect_error(doubletCluster(counts[,0], clusters[0]), "need at least three")
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

###################################################

set.seed(9900002)
test_that("doubletCells PC spawning works correctly", {
    sf <- runif(ncol(counts))
    y <- log2(t(t(counts)/sf)+1)
    centers <- rowMeans(y)
    SVD <- svd(t(y - centers), nv=20)

    set.seed(12345)
    sim.pcs <- scran:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=10000L, block=10000L)

    set.seed(12345)
    L <- sample(ncol(counts), 10000L, replace=TRUE)
    R <- sample(ncol(counts), 10000L, replace=TRUE)
    ref.x <- counts[,L] + counts[,R]
    ref.y <- log2(t(t(ref.x)/(sf[L] + sf[R]))+1)
    ref.pcs <- crossprod(ref.y - centers, SVD$v)

    expect_equal(sim.pcs, ref.pcs)

    # Works with multiple iterations.
    set.seed(23456)
    sim.pcs <- scran:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=25000L, block=10000L)

    set.seed(23456)
    ref1 <- scran:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=10000L, block=10000L)
    ref2 <- scran:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=10000L, block=10000L)
    ref3 <- scran:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=5000L, block=10000L)

    expect_equal(sim.pcs, rbind(ref1, ref2, ref3))
    expect_identical(dim(sim.pcs), c(25000L, ncol(SVD$v)))
})

set.seed(9900003)
test_that("size factor variations in doubletCells work correctly", {
    sf1 <- runif(ncol(counts))
    sf2 <- runif(ncol(counts))

    set.seed(12345)
    out <- doubletCells(counts)
    set.seed(12345)
    ref <- doubletCells(counts, size.factors.norm=scater::librarySizeFactors(counts))
    expect_equal(out, ref)

    set.seed(23456)
    out <- doubletCells(counts, size.factors.norm=sf1)
    set.seed(23456)
    ref <- doubletCells(t(t(counts)/sf1), size.factors.norm=rep(1, ncol(counts)), size.factors.content=1/sf1)
    expect_equal(out, ref)

    set.seed(34567)
    out <- doubletCells(counts, size.factors.norm=sf1, size.factors.content=sf2)
    set.seed(34567)
    ref <- doubletCells(t(t(counts)/sf2), size.factors.norm=sf1/sf2)
    expect_equal(out, ref)
})

set.seed(9900004)
test_that("high-level tests for doubletCells work correctly", {
    mu1 <- 2^rnorm(ngenes)
    mu2 <- 2^rnorm(ngenes)
    ncA <- 100
    ncB <- 100
    ncC <- 20

    counts.A <- matrix(mu1, ncol=ncA, nrow=ngenes)
    counts.B <- matrix(mu2, ncol=ncB, nrow=ngenes)
    counts.C <- matrix(mu1+mu2, ncol=ncC, nrow=ngenes)
    clusters <- rep(1:3, c(ncA, ncB, ncC))

    out <- doubletCells(cbind(counts.A, counts.B, counts.C))
    expect_true(min(out[clusters==3]) > max(out[clusters!=3]))

    # Now with differences in RNA content.
    counts.A <- matrix(mu1, ncol=ncA, nrow=ngenes)
    counts.B <- matrix(mu2, ncol=ncB, nrow=ngenes)
    counts.C <- matrix(mu1+2*mu2, ncol=ncC, nrow=ngenes)
    sf.spike <- 1/rep(1:3, c(ncA, ncB, ncC))
    
    X <- cbind(counts.A, counts.B, counts.C) 
    out <- doubletCells(X, size.factors.content=sf.spike)
    expect_true(min(out[clusters==3]) > max(out[clusters!=3]))
    expect_true(min(out[clusters==3]) > 2 * max(out[clusters!=3]))

    out <- doubletCells(X) # fails without size factor info; differences are basically negligible.
    expect_true(max(out[clusters==3]) < 2 * min(out[clusters!=3]))

    out <- scran:::.doublet_cells(X, force.match=TRUE, k=20) # recovers with forced matching.
    expect_true(min(out[clusters==3]) > max(out[clusters!=3]))
    expect_true(min(out[clusters==3]) > 2*max(out[clusters!=3]))
})

set.seed(9900005)
test_that("other settings for doubletCells work correctly", {
    # Subsetting behaves correctly.
    set.seed(1000)
    sim <- doubletCells(counts, subset.row=1:50)
    set.seed(1000)
    ref <- doubletCells(counts[1:50,])
    expect_identical(sim, ref)

    # Warnings raised if too many neighbors are requested.
    expect_warning(doubletCells(counts, k=1000), "'k' capped")

    # IRLBA works correctly.
    set.seed(2000)
    sim <- doubletCells(counts, d=5)
    set.seed(2000)
    ref <- doubletCells(counts, approximate=TRUE, irlba.args=list(tol=1e-12, work=50, maxit=20000), d=5)
    expect_true(median( abs(sim-ref)/(sim+ref+1e-6) ) < 0.01)

    # Responds correctly to blocking.
    set.seed(3000)
    ref <- doubletCells(counts)
    sim1 <- doubletCells(counts, block=1000)
    expect_equal(sim1, ref, tol=0.1)
    sim2 <- doubletCells(counts, niters=20000)
    expect_equal(sim2, ref, tol=0.1)
})

set.seed(9900006)
test_that("doubletCells works correctly for SCE objects", {
    sce <- SingleCellExperiment(list(counts=counts))

    set.seed(1000)
    ref <- doubletCells(counts)
    set.seed(1000)
    dbl <- doubletCells(sce)
    expect_identical(ref, dbl)

    # With a different assay.
    assay(sce, "whee") <- counts + rpois(length(counts), lambda=2)
    set.seed(1001)
    ref2 <- doubletCells(assay(sce, "whee"))
    set.seed(1001)
    dbl2 <- doubletCells(sce, assay.type="whee")
    expect_identical(ref2, dbl2)

    # With spike-ins that get used.
    isSpike(sce, "ERCC") <- sample(nrow(sce), 20)
    set.seed(1000)
    dbl3 <- doubletCells(sce, get.spikes=TRUE)
    expect_identical(ref, dbl3)

    # ... or ignored.
    set.seed(1002)
    dbl4 <- doubletCells(sce)
    set.seed(1002)
    ref3 <- doubletCells(counts(sce), subset.row=!isSpike(sce))
    expect_identical(ref3, dbl4)

    # With both spike-ins _and_ subset.row specified.
    rand.nospiked <- sample(which(!isSpike(sce)), 10)
    rand.spiked <- sample(which(isSpike(sce)), 10)

    set.seed(1003)
    expect_warning(dbl5 <- doubletCells(sce, subset.row=c(rand.spiked, rand.nospiked)), "greater than available")
    set.seed(1003)
    expect_warning(ref4 <- doubletCells(counts(sce), subset.row=rand.nospiked),  "greater than available")
    expect_identical(ref4, dbl5)
})


