# This checks the correlatePairs function.
# require(scran); require(testthat); source("setup.R"); source("test-correlate-pairs.R")

checkCorrelations <- function(out, exprs, null.dist) {
    ranked.exprs <- apply(exprs, 1, FUN=rank, ties.method="average")

    all.ranks <- seq_len(ncol(exprs))
    ranked.exprs <- ranked.exprs - mean(all.ranks)
    mean.sqdiff <- mean((all.ranks - mean(all.ranks))^2)
    ranked.exprs <- ranked.exprs / sqrt(mean.sqdiff)

    colnames(ranked.exprs) <- rownames(exprs)
    assembled.pval <- assembled.rho <- numeric(nrow(out))
    for (p in seq_along(assembled.rho)) { 
        assembled.rho[p] <- mean(ranked.exprs[,out$gene1[p]] * ranked.exprs[,out$gene2[p]])
        assembled.pval[p] <- min(sum(null.dist <= assembled.rho[p] + 1e-8), sum(null.dist >= assembled.rho[p] - 1e-8))
    }

    assembled.pval <- 2*(assembled.pval + 1)/(length(null.dist)+1)
    assembled.pval <- pmin(assembled.pval, 1)
    data.frame(rho=assembled.rho, pvalue=assembled.pval, FDR=p.adjust(assembled.pval, method="BH"))
}

set.seed(10000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

test_that("correlatePairs works with or without a pre-specified null distribution", {
    nulls <- sort(runif(1e5, -1, 1))

    out <- correlatePairs(X, null.dist=nulls)
    ref <- checkCorrelations(out, X, null.dist=nulls)
    
    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)
 
    set.seed(100)
    out <- correlatePairs(X, iters=1e5)
    set.seed(100)
    nulls <- correlateNull(Ncells, iter=1e5)
    ref <- checkCorrelations(out, X, null.dist=nulls)
    
    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)
})

set.seed(100001)
test_that("correlatePairs identifies limited pairs correctly", {
    X <- rbind(1:Ncells, 1:Ncells, as.numeric(rbind(1:(Ncells/2), Ncells - 1:(Ncells/2) + 1L)))
    out <- correlatePairs(X, null.dist=runif(1000))
    expect_identical(out$gene1, c(1L, 1L, 2L))
    expect_identical(out$gene2, c(2L, 3L, 3L))
    expect_identical(out$limited, c(TRUE, FALSE, FALSE))
})

set.seed(100002)
test_that("standard rho calculation works correctly", {
    nulls <- sort(runif(1e5, -1, 1))

    old <- correlatePairs(X, null.dist=nulls)
    out <- correlatePairs(X, null.dist=nulls, ties.method="average")
    expect_false(isTRUE(all.equal(old$rho, out$rho)))

    assembled.pval <- assembled.rho <- numeric(nrow(out))
    for (p in seq_along(assembled.rho)) { 
        assembled.rho[p] <- cor(X[out$gene1[p],], X[out$gene2[p],], method="spearman")
        assembled.pval[p] <- min(sum(nulls <= assembled.rho[p] + 1e-8), sum(nulls >= assembled.rho[p] - 1e-8))
    }

    assembled.pval <- 2*(assembled.pval + 1)/(length(nulls)+1)
    assembled.pval <- pmin(assembled.pval, 1)

    expect_equal(out$rho, assembled.rho)
    expect_equal(out$p.value, assembled.pval)

    # No difference if there aren't ties in the first place.
    Y <- matrix(rnorm(1000), ncol=100)
    old <- correlatePairs(Y, null.dist=nulls)
    out <- correlatePairs(Y, null.dist=nulls, ties.method="average")
    expect_equal(old, out)
})

####################################################################################################

test_that("correlatePairs works with a design matrix", {
    grouping <- gl(2, 50)
    design <- model.matrix(~grouping)

    # Test the residual calculator.
    QR <- qr(design, LAPACK=TRUE)
    ref.resid <- t(lm.fit(y=t(X), x=design)$residuals)
    out.resid <- scran:::get_residuals(X, QR$qr, QR$qraux, seq_len(nrow(X))-1L, NA_real_) 
    expect_equal(unname(ref.resid), out.resid)
    
    subset.chosen <- sample(nrow(X), 10)
    out.resid.sub <- scran:::get_residuals(X, QR$qr, QR$qraux, subset.chosen-1L, NA_real_) 
    expect_equal(unname(ref.resid[subset.chosen,]), out.resid.sub)

    # Manual calculation (note; using out.resid as ties are slightly different with ref.resid).
    set.seed(200)
    out <- correlatePairs(X, design=design, iters=1e4)

    set.seed(200)
    nulls <- correlateNull(design=design, iter=1e4)
    rownames(out.resid) <- rownames(X)
    ref <- checkCorrelations(out, out.resid, null.dist=nulls)

    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)
})

test_that("correlatePairs works with blocking", {
    grouping <- sample(LETTERS[1:3], Ncells, replace=TRUE)

    # Need to use a normal matrix to avoid ties; bplapply mucks up 
    # the seed for tie breaking, even with a single core.
    X <- matrix(rnorm(Ncells*Ngenes), ncol=Ncells)

    # Without weighting.
    set.seed(200)
    out <- correlatePairs(X, block=grouping, iters=1e4)

    set.seed(200)
    nulls <- correlateNull(block=grouping, iter=1e4)
    collected.rho <- 0L
    for (group in split(seq_along(grouping), grouping)) { 
        ref <- checkCorrelations(out, X[,group], null.dist=nulls)
        collected.rho <- collected.rho + ref$rho/length(unique(grouping))
    }
    expect_equal(out$rho, collected.rho)

    # Obtaining a p-value.
    collected.p <- numeric(length(collected.rho)) 
    for (x in seq_along(collected.rho)) { 
        collected.p[x] <- min(sum(nulls <= collected.rho[x] + 1e-8), sum(nulls >= collected.rho[x] - 1e-8))
    }
    collected.p <- 2*(collected.p + 1)/(length(nulls)+1)
    collected.p <- pmin(collected.p, 1)
    expect_equal(out$p.value, collected.p)

    # With weighting.
    set.seed(200)
    out <- correlatePairs(X, block=grouping, equiweight=FALSE, iters=1e4)

    set.seed(200)
    nulls <- correlateNull(block=grouping, equiweight=FALSE, iter=1e4)
    collected.rho <- 0L
    for (group in split(seq_along(grouping), grouping)) { 
        ref <- checkCorrelations(out, X[,group], null.dist=nulls)
        collected.rho <- collected.rho + length(group)/length(grouping) * ref$rho
    }
    expect_equal(out$rho, collected.rho)
    
    # Obtaining a p-value.
    collected.p <- numeric(length(collected.rho)) 
    for (x in seq_along(collected.rho)) { 
        collected.p[x] <- min(sum(nulls <= collected.rho[x] + 1e-8), sum(nulls >= collected.rho[x] - 1e-8))
    }
    collected.p <- 2*(collected.p + 1)/(length(nulls)+1)
    collected.p <- pmin(collected.p, 1)
    expect_equal(out$p.value, collected.p)
})

####################################################################################################
# Checking that it works with different 'pairing' values.
# (again, using a normal matrix to avoid ties, for simplicity).

set.seed(10002)
Ngenes <- 20
Ncells <- 100
X <- matrix(rnorm(Ngenes*Ncells), nrow=Ngenes)
rownames(X) <- paste0("X", seq_len(Ngenes))

nulls <- correlateNull(ncells=ncol(X), iter=1e4)
ref <- correlatePairs(X, null.dist=nulls)

set.seed(100021)
test_that("correlatePairs works with subset.row values", {
    # Vanilla subsetting of genes; taking the correlations between all pairs of subsetted genes. 
    subgenes <- 1:10 
    subbed <- correlatePairs(X, nulls, subset.row=subgenes)
    expected <- correlatePairs(X[subgenes,], nulls)
    expect_equal(subbed, expected)

    # More rigorous check relative to the reference.
    expect_true(!is.unsorted(subbed$p.value))
    subref <- ref[ref$gene1 %in% rownames(X)[subgenes] & ref$gene2 %in% rownames(X)[subgenes],]
    subref$FDR <- p.adjust(subref$p.value, method="BH")
    rownames(subref) <- NULL
    expect_equal(subref, subbed)
})

set.seed(100022)
test_that("correlatePairs works with list-based pairing", {
    ref2 <- ref
    ref2$gene1 <- ref$gene2
    ref2$gene2 <- ref$gene1
    ref.com <- rbind(ref, ref2)
    ref.names <- paste0(ref.com$gene1, ".", ref.com$gene2)

    for (attempt in 1:3) {
        if (attempt==1L) {
            pairs <- list(1:10, 5:15)
        } else if (attempt==2L) {
            pairs <- list(1:10, 11:20)
        } else {
            pairs <- list(20:10, 15:11)
        }
    
        lsubbed <- correlatePairs(X, nulls, pairings=pairs)
        expect_true(!is.unsorted(lsubbed$p.value))
        expect_true(all(lsubbed$gene1!=lsubbed$gene2) && all(lsubbed$gene1 %in% rownames(X)[pairs[[1]]])
                    && all(lsubbed$gene2 %in% rownames(X)[pairs[[2]]]))
        overlap <- length(intersect(pairs[[1]], pairs[[2]]))
        expect_true(nrow(lsubbed) == length(pairs[[1]]) * length(pairs[[2]]) - overlap)
    
        lsub.names <- paste0(lsubbed$gene1, ".", lsubbed$gene2)
        expect_true(!any(duplicated(lsub.names)))
        lsub.ref <- ref.com[match(lsub.names, ref.names),]
        lsub.ref$FDR <- p.adjust(lsub.ref$p.value, method="BH")
        rownames(lsub.ref) <- NULL
        expect_equal(lsub.ref, lsubbed)

        # Checking for proper interaction with subset.row
        keep <- rep(c(TRUE, FALSE), length.out=nrow(X))
        lsub2 <- correlatePairs(X, nulls, pairings=pairs, subset.row=keep)
        repairs <- list(intersect(pairs[[1]], which(keep)), intersect(pairs[[2]], which(keep)))
        lexpected <- correlatePairs(X, nulls, pairings=repairs)
        expect_equal(lexpected, lsub2)
    }

    expect_error(correlatePairs(X, nulls, pairings=list()), "'pairings' as a list should have length 2")
    expect_error(correlatePairs(X, nulls, pairings=list(1,2,3)), "'pairings' as a list should have length 2")
})

set.seed(100023)
test_that("correlatePairs with pairs matrix works as expected", {
    ref2 <- ref
    ref2$gene1 <- ref$gene2
    ref2$gene2 <- ref$gene1
    ref.com <- rbind(ref, ref2)
    ref.names <- paste0(ref.com$gene1, ".", ref.com$gene2)

    for (attempt in 1:3) { 
        if (attempt==1L) {
            pairs <- cbind(1:10, 2:11)
        } else if (attempt==2L) {
            pairs <- cbind(20:6, 5:19)
        } else {
            choices <- which(!diag(nrow(X)), arr.ind=TRUE)
            pairs <- choices[sample(nrow(choices), 20),]
        }

        msubbed <- correlatePairs(X, nulls, pairings=pairs)
        expect_identical(rownames(X)[pairs[,1]], msubbed$gene1)
        expect_identical(rownames(X)[pairs[,2]], msubbed$gene2)
    
        msub.names <- paste0(msubbed$gene1, ".", msubbed$gene2)
        msub.ref <- ref.com[match(msub.names, ref.names),]
        msub.ref$FDR <- p.adjust(msub.ref$p.value, method="BH")
        rownames(msub.ref) <- NULL
        expect_equal(msub.ref, msubbed)
    
        pairs2 <- rownames(X)[pairs]
        dim(pairs2) <- dim(pairs)
        msubbed2 <- correlatePairs(X, nulls, pairings=pairs2)
        expect_equal(msubbed, msubbed2)
    
        # Checking for proper interaction with subset.row
        keep <- 1:10
        lsub2 <- correlatePairs(X, nulls, pairings=pairs, subset.row=keep)
        repairs <- pairs[pairs[,1] %in% keep & pairs[,2] %in% keep,,drop=FALSE]
        lexpected <- correlatePairs(X, nulls, pairings=repairs)
        expect_equal(lexpected, lsub2)
    }

    expect_error(correlatePairs(X, nulls, pairings=matrix(0, 0, 0)), 
        "'pairings' should be a numeric/character matrix with 2 columns")
    expect_error(correlatePairs(X, nulls, pairings=cbind(1,2,3)), 
        "'pairings' should be a numeric/character matrix with 2 columns")
    expect_error(correlatePairs(X, nulls, pairings=cbind(TRUE, FALSE)), 
        "'pairings' should be a numeric/character matrix with 2 columns")
})

####################################################################################################

set.seed(10003)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))
nulls <- sort(runif(1e6, -1, 1))

test_that("correlatePairs works correctly with SingleCellExperiment objects", {
    set.seed(100)
    ref <- correlatePairs(X, null.dist=nulls)
    set.seed(100)
    X2 <- SingleCellExperiment(list(logcounts=X))
    metadata(X2)$log.exprs.offset <- 1
    out <- correlatePairs(X2, null.dist=nulls)
    expect_equal(out, ref)
    
    # With spikes.
    isSpike(X2, "MySpike") <- rbinom(Ngenes, 1, 0.6)==0L
    set.seed(100)
    ref <- correlatePairs(exprs(X2)[!isSpike(X2),,drop=FALSE], null.dist=nulls)
    set.seed(100)
    out <- correlatePairs(X2, null.dist=nulls)
    expect_equal(out, ref)
})

test_that("correlatePairs fails properly upon silly inputs", {
    expect_error(correlatePairs(X[0,], nulls), "need at least two genes to compute correlations")
    expect_true(all(is.nan(correlatePairs(X[,0], nulls)$rho)))
    expect_true(all(is.na(correlatePairs(X, nulls, design=diag(ncol(X)))$rho)))
    expect_true(all(is.nan(correlatePairs(X, nulls, block=1:ncol(X))$rho)))

    expect_warning(correlatePairs(X, iters=1), "lower bound on p-values at a FDR of 0.05, increase 'iter'")
    expect_warning(out <- correlatePairs(X, numeric(0)), "lower bound on p-values at a FDR of 0.05, increase 'iter'")
    expect_equal(out$p.value, rep(1, nrow(out)))
})
