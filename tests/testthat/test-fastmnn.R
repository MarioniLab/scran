# This tests the functions related to fastMNN.
# library(scran); library(testthat); source("test-fastmnn.R")

set.seed(1200001)
test_that("averaging correction vectors works as expected", {
    test1 <- matrix(rnorm(1000), ncol=10)
    test2 <- matrix(rnorm(2000), ncol=10)
    mnn1 <- sample(nrow(test1), 250, replace=TRUE)
    mnn2 <- sample(nrow(test1), 250, replace=TRUE)

    # Slow reference calculation.
    correct <- test1[mnn1,] - test2[mnn2,]
    by.mnn <- split(seq_along(mnn2), mnn2)
    collected <- vector("list", length(by.mnn))
    for (idx in seq_along(by.mnn)) {
        collected[[idx]] <- colMeans(correct[by.mnn[[idx]],,drop=FALSE])
    }
    ref <- do.call(rbind, collected)
    rownames(ref) <- names(by.mnn)

    # Comparing the implementation in scran.
    out <- scran:::.average_correction(test1, mnn1, test2, mnn2)  
    expect_equal(out$averaged, ref)
    expect_identical(out$second, sort(unique(mnn2)))
    expect_identical(as.integer(rownames(out$averaged)), out$second)
})

set.seed(1200002)
test_that("centering along a batch vector works correctly", {
    test <- matrix(rnorm(1000), ncol=10)
    batch <- rnorm(10) 
    centered <- scran:::.center_along_batch_vector(test, batch)
    new.locations <- centered %*% batch
    expect_true(mad(new.locations) < 1e-8)
})

set.seed(1200003)
test_that("tricube weighting works correctly", {
    test <- matrix(rnorm(1000), ncol=10)
    correction <- matrix(rnorm(500), ncol=10)
    involved <- sample(nrow(test), nrow(correction))

    # Setting up a reference function for comparison, operating truly row-by-row.
    FUN <- function(current, corvec, in.mnn, k=20, ndist=3) {
        cur.uniq <- current[in.mnn,,drop=FALSE]
        safe.k <- min(k, nrow(cur.uniq))
        closest <- kmknn::queryKNN(query=current, X=cur.uniq, k=safe.k)
        middle.k <- ceiling(safe.k/2L)

        for (x in seq_len(nrow(current))) {
            all.dists <- closest$distance[x,]
            all.index <- closest$index[x,]

            middist <- sort(all.dists)[middle.k] 
            weights <- (1 - pmin(1, all.dists/(middist*ndist))^3)^3
            weights <- weights/sum(weights)

            curcor <- colSums(corvec[all.index,] * weights) 
            current[x,] <- current[x,] + curcor
        }

        return(current)
    }

    out <- scran:::.tricube_weighted_correction(test, correction, involved, k=20, ndist=3)
    ref <- FUN(test, correction, involved, k=20, ndist=3)
    expect_equal(ref, out)

    out <- scran:::.tricube_weighted_correction(test, correction, involved, k=11, ndist=3)
    ref <- FUN(test, correction, involved, k=11, ndist=3)
    expect_equal(ref, out)

    out <- scran:::.tricube_weighted_correction(test, correction, involved, k=11, ndist=1)
    ref <- FUN(test, correction, involved, k=11, ndist=1)
    expect_equal(ref, out)
})

CHECK_PAIRINGS <- function(origin, pairings) {
    expect_identical(unique(origin$batch[pairings[[1]]$first]), origin$batch[1])

    for (p in pairings) {
        expect_true(all(p$first < p$second))

        sbatch <- origin$batch[p$second]
        expect_identical(length(unique(sbatch)), 1L)

        fbatch <- origin$batch[p$first]
        expect_true(!any(fbatch %in% sbatch))
    }

    return(NULL)
}

set.seed(1200004)
test_that("fastMNN works as expected for two batches", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, d=50) # corrected values
    expect_identical(dim(out$corrected), c(ncol(B1) + ncol(B2), 50L))
    expect_identical(out$origin$batch, rep(1:2, c(ncol(B1), ncol(B2))))
    expect_identical(out$origin$cell, c(seq_len(ncol(B1)), seq_len(ncol(B2))))
    CHECK_PAIRINGS(out$origin, out$pairs)
    
    # Dimension choice behaves correctly.
    out.10 <- fastMNN(B1, B2, d=10) 
    expect_identical(ncol(out.10$corrected), 10L)
    CHECK_PAIRINGS(out.10$origin, out.10$pairs)

    # Handles names correctly.
    out.n <- fastMNN(X=B1, Y=B2, d=50) 
    expect_identical(out$corrected, out.n$corrected)
    expect_identical(out.n$origin$batch, rep(c("X", "Y"), c(ncol(B1), ncol(B2))))
    expect_identical(out.n$origin$cell, c(seq_len(ncol(B1)), seq_len(ncol(B2))))
    CHECK_PAIRINGS(out.n$origin, out.n$pairs)

    # Behaves if we turn off cosine-norm.
    nB1 <- t(t(B1)/ sqrt(colSums(B1^2)))
    nB2 <- t(t(B2)/ sqrt(colSums(B2^2)))
    out.ncos <- fastMNN(nB1, nB2, cos.norm=FALSE, d=50) 
    expect_equal(out.ncos, out) 

    # Subset.row behaves correctly.
    i <- sample(nrow(B1), 50)
    ref <- fastMNN(X=B1[i,], Y=B2[i,], d=50)
    out.s <- fastMNN(X=B1, Y=B2, d=50, subset.row=i)
    expect_identical(out.s, ref)

    # Behaves if we only use PCs.
    pcs <- multiBatchPCA(B1, B2, d=10, approximate=FALSE)
    out.pre <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
    out.norm <- fastMNN(B1, B2, d=10, cos.norm=FALSE, approximate=FALSE)
    expect_equal(out.pre, out.norm)
})

set.seed(1200005)
test_that("fastMNN works as expected for three batches, with auto-ordering", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, B3, d=50) # corrected values
    expect_identical(dim(out$corrected), c(ncol(B1) + ncol(B2) + ncol(B3), 50L))
    expect_identical(out$origin$batch, rep(1:3, c(ncol(B1), ncol(B2), ncol(B3))))
    expect_identical(out$origin$cell, c(seq_len(ncol(B1)), seq_len(ncol(B2)), seq_len(ncol(B3))))
    CHECK_PAIRINGS(out$origin, out$pairs)

    # Testing the auto-ordering algorithms. 
    out.auto <- fastMNN(B1, B2, B3, d=50, auto.order=TRUE) 
    expect_identical(dim(out$corrected), dim(out.auto$corrected))
    expect_identical(out.auto$origin$batch, rep(c(2L, 1L, 3L), c(ncol(B2), ncol(B1), ncol(B3)))) # 3 should be last, with the fewest cells => fewest MNNs.
    expect_identical(out.auto$origin$cell, c(seq_len(ncol(B2)), seq_len(ncol(B1)), seq_len(ncol(B3))))
    CHECK_PAIRINGS(out.auto$origin, out.auto$pairs)

    # Testing the internal auto-ordering functions.
    fmerge <- scran:::.define_first_merge(list(t(B1), t(B2), t(B3)), k=20)   
    expect_identical(fmerge$first, 2L) # 2 is arbitrarily 'first', and 1 is arbitrarily 'second'; but 3 should never show up.
    expect_identical(fmerge$second, 1L)
    expect_identical(fmerge$pairs, scran:::find.mutual.nn(t(B2), t(B1), k1=20, k2=20, BPPARAM=SerialParam()))

    expect_identical(kmknn::findKNN(precomputed=fmerge$precomputed[[1]], k=5), kmknn::findKNN(t(B1), k=5)) # checking that the precomputations are correct.
    expect_identical(kmknn::findKNN(precomputed=fmerge$precomputed[[2]], k=10), kmknn::findKNN(t(B2), k=10))
    expect_identical(kmknn::findKNN(precomputed=fmerge$precomputed[[3]], k=15), kmknn::findKNN(t(B3), k=15))

    nmerge <- scran:::.define_next_merge(t(B1), list(t(B1), t(B2), t(B3)), processed=1L, precomputed=fmerge$precomputed, k=20)   
    expect_identical(nmerge$other, 2L) # 1 should have more MNNs with 2 than with 3.
    expect_identical(nmerge$pairs, scran:::find.mutual.nn(t(B1), t(B2), k1=20, k2=20, BPPARAM=SerialParam()))
})

set.seed(1200006)
test_that("fastMNN fails on silly inputs", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    # Throws errors properly with no genes or no cells.
    expect_error(fastMNN(), "at least two batches")
    expect_error(fastMNN(B1), "at least two batches")
    expect_error(fastMNN(B1[0,], B2[0,]), "too large")
    expect_error(fastMNN(B1[,0], B2[,0]), "too large")

    # Throws errors upon row checks.
    expect_error(fastMNN(B1[1:10,], B2), "number of rows is not the same")
    xB1 <- B1
    xB2 <- B2
    rownames(xB1) <- sample(nrow(B1))
    rownames(xB2) <- sample(nrow(B2))
    expect_error(fastMNN(xB1, xB2), "row names are not the same")
})
