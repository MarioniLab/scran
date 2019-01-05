# This checks the correlateNull function.
# require(scran); require(testthat); source("setup.R"); source("test-correlate-pairs.R")

.tolerant_rank <- function(y, tol=1e-6) {
    if (!length(y)) { return(integer(0)) }
    o <- order(y)
    rle.out <- rle(y[o])
    okay <- c(TRUE, diff(rle.out$values) > tol)
    to.use <- cumsum(okay)
    rle.out$values <- rle.out$values[okay][to.use]
    y[o] <- inverse.rle(rle.out)
    rank(y, ties.method="random")
}

set.seed(1000)
test_that("error tolerant ranking is working correctly", {
    whee <- runif(100, -1e-16, 1e-16)
    set.seed(100)
    r <- .tolerant_rank(whee)
    set.seed(100)
    r2 <- rank(integer(100), ties.method="random")
    set.seed(100)
    r3 <- .Call(scran:::cxx_get_untied_ranks, rbind(whee), 0L, seq_along(whee)-1L, 1e-6)
    
    expect_identical(r, r2)
    expect_identical(r, r3[,1])
    
    set.seed(200)
    extra <- sample(10, 100, replace=TRUE)
    set.seed(100)
    r <- .tolerant_rank(whee + extra)
    set.seed(100)
    r2 <- rank(extra, ties.method="random")
    set.seed(100)
    r3 <- .Call(scran:::cxx_get_untied_ranks, rbind(whee + extra), 0L, seq_along(extra)-1L, 1e-6)
    set.seed(100)
    r4 <- .Call(scran:::cxx_get_untied_ranks, rbind(extra), 0L, seq_along(extra)-1L, 1e-6)
    
    expect_identical(r, r2)
    expect_identical(r, r3[,1])
    expect_identical(r, r4[,1])
})

####################################################################################################

checkCorrelations <- function(out, exprs, null.dist) {
    ranked.exprs <- apply(exprs, 1, FUN=.tolerant_rank)
    colnames(ranked.exprs) <- rownames(exprs)
    assembled.pval <- assembled.rho <- numeric(nrow(out))
    for (p in seq_along(assembled.rho)) { 
        assembled.rho[p] <- cor(ranked.exprs[,out$gene1[p]], ranked.exprs[,out$gene2[p]], method="spearman")        
        assembled.pval[p] <- min(sum(null.dist <= assembled.rho[p] + 1e-8), sum(null.dist >= assembled.rho[p] - 1e-8))
    }
    assembled.pval <- 2*(assembled.pval + 1)/(length(null.dist)+1)
    assembled.pval <- pmin(assembled.pval, 1)
    return(data.frame(rho=assembled.rho, pvalue=assembled.pval, FDR=p.adjust(assembled.pval, method="BH")))
}

set.seed(10000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

test_that("correlatePairs works with or without a pre-specified null distribution", {
    nulls <- sort(runif(1e5, -1, 1))

    set.seed(100)
    out <- correlatePairs(X, null.dist=nulls, tie.iters=1)
    set.seed(100)
    ref <- checkCorrelations(out, X, null.dist=nulls)
    
    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)
 
    set.seed(100)
    out <- correlatePairs(X, null.iters=1e5, tie.iters=1)
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
test_that("tie cracking works correctly", {
    nulls <- sort(runif(1e5, -1, 1))

    set.seed(100)
    out1 <- correlatePairs(X, null.dist=nulls, tie.iters=1)
    out2 <- correlatePairs(X, null.dist=nulls, tie.iters=1)
    out3 <- correlatePairs(X, null.dist=nulls, tie.iters=1)

    out1 <- out1[order(out1$gene1, out1$gene2),]
    out2 <- out2[order(out2$gene1, out2$gene2),]
    out3 <- out3[order(out3$gene1, out3$gene2),]

    # Checking that the random seed does change.
    expect_false(isTRUE(all.equal(out1$rho, out2$rho)))
    expect_false(isTRUE(all.equal(out1$rho, out3$rho)))

    # Checking that combined stats are calculated correctly.
    set.seed(100)
    out <- correlatePairs(X, null.dist=nulls, tie.iters=3)
    out <- out[order(out$gene1, out$gene2),]

    expect_equal(out$rho, (out1$rho + out2$rho + out3$rho)/3)
    expect_equal(out$p.value, combinePValues(out1$p.value, out2$p.value, out3$p.value, method="simes"))
    expect_identical(out$limited, out1$limited | out2$limited | out3$limited)
})

####################################################################################################

grouping <- gl(2, 50)
design <- model.matrix(~grouping)

test_that("correlatePairs works with a design matrix", {
    set.seed(200)
    nulls <- correlateNull(design=design, iter=1e4)
    expect_warning(correlatePairs(X[1:5,], design=NULL, null.dist=nulls), "'design' is not the same as that used to generate")
    
    set.seed(100) # Need because of random ranking.
    out <- correlatePairs(X, design=design, null.dist=nulls, lower.bound=NA, tie.iters=1)
    fit <- lm.fit(x=design, y=t(X))
    exprs <- t(fit$residual)
    set.seed(100)
    ref <- checkCorrelations(out, exprs, null.dist=nulls)
    
    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)
    
    # Deeper test of the residual calculator.
    QR <- qr(design, LAPACK=TRUE)
    ref.resid <- t(lm.fit(y=t(X), x=design)$residuals)
    out.resid <- .Call(scran:::cxx_get_residuals, X, QR$qr, QR$qraux, seq_len(nrow(X))-1L, NA_real_) 
    expect_equal(unname(ref.resid), out.resid)
    
    subset.chosen <- sample(nrow(X), 10)
    out.resid <- .Call(scran:::cxx_get_residuals, X, QR$qr, QR$qraux, subset.chosen-1L, NA_real_) 
    expect_equal(unname(ref.resid[subset.chosen,]), out.resid)
})

test_that("correlatePairs works with lower bounds on the residuals", {
    set.seed(100021)
    X[] <- log(matrix(rpois(Ngenes*Ncells, lambda=1), nrow=Ngenes) + 1) # more zeroes with lambda=1
    expect_true(sum(X==0)>0) # observations at the lower bound.
    
    # Checking proper calculation of the bounded residuals.
    ref.resid <- t(lm.fit(y=t(X), x=design)$residuals)
    alt.resid <- ref.resid
    is.smaller <- X <= 0 # with bounds
    smallest.per.row <- apply(alt.resid, 1, min) - 1
    alt.resid[is.smaller] <- smallest.per.row[arrayInd(which(is.smaller), dim(is.smaller))[,1]]

    QR <- qr(design, LAPACK=TRUE)
    out.resid <- .Call(scran:::cxx_get_residuals, X, QR$qr, QR$qraux, seq_len(nrow(X))-1L, 0)
    expect_equal(unname(alt.resid), out.resid)
   
    # Checking for correct calculation of statistics for bounded residuals.
    set.seed(200)
    nulls <- correlateNull(design=design, iter=1e4)
 
    set.seed(100) 
    out <- correlatePairs(X, design=design, null.dist=nulls, lower.bound=0, tie.iters=1)
    set.seed(100) 
    ref <- checkCorrelations(out, alt.resid, null.dist=nulls)
    
    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)
})

test_that("correlatePairs works without simulated residuals for one-way layouts", {
    # Repeating without simulation (need to use a normal matrix to avoid ties;
    # bplapply mucks up the seed for tie breaking, even with a single core).
    set.seed(200)
    X[] <- rnorm(length(X))
    nulls <- correlateNull(block=grouping, iter=1e4)
    out <- correlatePairs(X, block=grouping, null.dist=nulls, tie.iters=1)
    expect_warning(correlatePairs(X, block=NULL, null.dist=nulls), "'block' is not the same")
 
    # Calculating the weighted average of correlations. 
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

    expect_error(correlatePairs(X, nulls, pairings=matrix(0, 0, 0)), "'pairings' should be a numeric/character matrix with 2 columns")
    expect_error(correlatePairs(X, nulls, pairings=cbind(1,2,3)), "'pairings' should be a numeric/character matrix with 2 columns")
    expect_error(correlatePairs(X, nulls, pairings=cbind(TRUE, FALSE)), "'pairings' should be a numeric/character matrix with 2 columns")
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
    
    # With per.gene=TRUE.
    set.seed(100)
    ref <- correlatePairs(exprs(X2), null.dist=nulls, per.gene=TRUE)
    set.seed(100)
    out <- correlatePairs(X2, null.dist=nulls, per.gene=TRUE)
    expect_equal(out, ref)
    
    # With spikes.
    isSpike(X2, "MySpike") <- rbinom(Ngenes, 1, 0.6)==0L
    set.seed(100)
    ref <- correlatePairs(exprs(X2)[!isSpike(X2),,drop=FALSE], null.dist=nulls)
    set.seed(100)
    out <- correlatePairs(X2, null.dist=nulls)
    expect_equal(out, ref)
})

# Checking nonsense inputs.

test_that("correlatePairs fails properly upon silly inputs", {
    expect_error(correlatePairs(X[0,], nulls), "need at least two genes to compute correlations")
    expect_error(correlatePairs(X[,0], nulls), "number of cells should be greater than or equal to 2")
    expect_warning(correlatePairs(X, iters=1), "lower bound on p-values at a FDR of 0.05, increase 'iter'")
    expect_warning(out <- correlatePairs(X, numeric(0)), "lower bound on p-values at a FDR of 0.05, increase 'iter'")
    expect_equal(out$p.value, rep(1, nrow(out)))
})
