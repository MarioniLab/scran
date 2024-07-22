# This checks the correlatePairs function.
# require(scran); require(testthat); source("setup.R"); source("test-correlate-pairs.R")

####################################################################################################
# This checks the basics of the correlation calculator.

test_that("rhoToPValue works as expected", {
    # Without ties.
    x1 <- rnorm(20)
    x2 <- rnorm(20)

    y <- cor.test(x1, x2, method="spearman", alternative="less", exact=FALSE)
    out <- rhoToPValue(cor(x1, x2, method="spearman"), n=length(x1), positive=FALSE)
    expect_equal(y$p.value, out)

    y <- cor.test(x1, x2, method="spearman", alternative="greater", exact=FALSE)
    out <- rhoToPValue(cor(x1, x2, method="spearman"), n=length(x1), positive=TRUE)
    expect_equal(y$p.value, out)

    y <- cor.test(x1, x2, method="spearman", exact=FALSE)
    out2 <- rhoToPValue(cor(x1, x2, method="spearman"), n=length(x1))
    expect_equal(y$p.value, pmin(out2[[1]], out2[[2]])*2)

    # Works with ties.
    x1 <- rpois(20, lambda=3)
    x2 <- rpois(20, lambda=3)

    y <- cor.test(x1, x2, method="spearman", alternative="less", exact=FALSE)
    out <- rhoToPValue(cor(x1, x2, method="spearman"), n=length(x1), positive=FALSE)
    expect_equal(y$p.value, out)

    y <- cor.test(x1, x2, method="spearman", alternative="greater", exact=FALSE)
    out <- rhoToPValue(cor(x1, x2, method="spearman"), n=length(x1), positive=TRUE)
    expect_equal(y$p.value, out)

    y <- cor.test(x1, x2, method="spearman", exact=FALSE)
    out2 <- rhoToPValue(cor(x1, x2, method="spearman"), n=length(x1))
    expect_equal(y$p.value, pmin(out2[[1]], out2[[2]])*2)

    # Vectorized properly.
    r <- runif(20, -1, 1)
    out <- rhoToPValue(r, n=length(x1))
    ref <- lapply(r, FUN=rhoToPValue, n=length(x1))
    ref <- do.call(mapply, c(list(c), ref, list(SIMPLIFY=FALSE)))
    expect_equal(out, ref)
})

test_that("basic correlation works correctly", {
    x <- matrix(rnorm(100), ncol=20)
    output <- scran:::.calculate_rho(x)
    ref <- vapply(seq_len(nrow(x)), function(i) {
        y <- x[i,]
        half <- length(y)/2
        cor(head(y, half), tail(y, half), method="spearman")
    }, 0)
    expect_equal(output, ref)

    # Handles ties like a champ.
    x <- matrix(rpois(100, 1), ncol=20)
    output <- scran:::.calculate_rho(x)
    ref <- vapply(seq_len(nrow(x)), function(i) {
        y <- x[i,]
        half <- length(y)/2
        cor(head(y, half), tail(y, half), method="spearman")
    }, 0)
    expect_equal(output, ref)

    x <- matrix(0, nrow=5, ncol=20)
    output <- scran:::.calculate_rho(x)
    expect_true(all(is.na(output)))
})

####################################################################################################

checkCorrelations <- function(out, exprs) {
    assembled.rho <- numeric(nrow(out))
    for (p in seq_along(assembled.rho)) { 
        assembled.rho[p] <- cor(exprs[out$gene1[p],], exprs[out$gene2[p],], method="spearman")
    }
    p.out <- rhoToPValue(assembled.rho, ncol(exprs))
    assembled.pval <- pmin(p.out$positive, p.out$negative) * 2
    data.frame(rho=assembled.rho, pvalue=assembled.pval, FDR=p.adjust(assembled.pval, method="BH"))
}

set.seed(10000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

test_that("correlatePairs works correctly in basic mode", {
    out <- correlatePairs(X)
    ref <- checkCorrelations(out, X)

    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)

    # Doesn't shit itself with all-zeroes.
    empties <- matrix(0, nrow=5, ncol=20)
    out <- correlatePairs(empties)
    expect_equal(nrow(out), choose(nrow(empties), 2L))
    expect_true(all(is.na(out$p.value)))
})

test_that("correlatePairs works with a design matrix", {
    grouping <- gl(2, 50)
    design <- model.matrix(~grouping)
    X <- matrix(rnorm(Ngenes*Ncells), nrow=Ngenes) # avoid problems with ties.

    out <- correlatePairs(X, design=design)
    resids <- t(lm.fit(y=t(X), x=design)$residuals)
    ref <- checkCorrelations(out, resids)

    expect_equal(out$rho, ref$rho)
    expect_equal(out$p.value, ref$pvalue)
    expect_equal(out$FDR, ref$FDR)

    # Works with subsetting.
    sub <- correlatePairs(X, design=design, subset.row=1:10)
    keep <- out$gene1 %in% 1:10 & out$gene2 %in% 1:10
    sub$FDR <- out$FDR <- NULL
    expect_identical(sub, out[keep,])
})

test_that("correlatePairs works with blocking", {
    grouping <- sample(LETTERS[1:3], Ncells, replace=TRUE)

    # Need to use a normal matrix to avoid ties; bplapply mucks up 
    # the seed for tie breaking, even with a single core.
    X <- matrix(rnorm(Ncells*Ngenes), ncol=Ncells)

    # Without weighting.
    out <- correlatePairs(X, block=grouping)

    by.group <- split(seq_along(grouping), grouping)
    all.rho <- vector("list", length(by.group))
    for (g in seq_along(by.group)) { 
        group <- by.group[[g]]
        ref <- checkCorrelations(out, X[,group])
        all.rho[[g]] <- ref$rho
    }
    expect_equal(out$rho, Reduce("+", all.rho)/length(by.group))

    # With weighting.
    out2 <- correlatePairs(X, block=grouping, equiweight=FALSE, pairings=as.matrix(out[,1:2]))
    scaled.rho <- mapply("*", all.rho, lengths(by.group), SIMPLIFY=FALSE)
    expect_equal(out2$rho, Reduce("+", scaled.rho)/ncol(X))

    # Handles failed blocks.
    failed <- grouping=="A"
    X[,failed] <- 0
    out <- correlatePairs(X, block=grouping)
    ref <- correlatePairs(X[,!failed], block=grouping[!failed])
    expect_equal(out, ref)
})

####################################################################################################
# Checking that it works with different 'pairing' values.
# (again, using a normal matrix to avoid ties, for simplicity).

set.seed(10002)
Ngenes <- 20
Ncells <- 100
X <- matrix(rnorm(Ngenes*Ncells), nrow=Ngenes)
rownames(X) <- paste0("X", seq_len(Ngenes))
ref <- correlatePairs(X)

set.seed(100021)
test_that("correlatePairs works with subset.row values", {
    # Vanilla subsetting of genes; taking the correlations between all pairs of subsetted genes. 
    subgenes <- 1:10 
    subbed <- correlatePairs(X, subset.row=subgenes)
    expected <- correlatePairs(X[subgenes,])
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
    
        lsubbed <- correlatePairs(X, pairings=pairs)
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
        lsub2 <- correlatePairs(X, pairings=pairs, subset.row=keep)
        repairs <- list(intersect(pairs[[1]], which(keep)), intersect(pairs[[2]], which(keep)))
        lexpected <- correlatePairs(X, pairings=repairs)
        expect_equal(lexpected, lsub2)
    }

    expect_error(correlatePairs(X, pairings=list()), "'pairings' as a list should have length 2")
    expect_error(correlatePairs(X, pairings=list(1,2,3)), "'pairings' as a list should have length 2")
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

        msubbed <- correlatePairs(X, pairings=pairs)
        expect_identical(rownames(X)[pairs[,1]], msubbed$gene1)
        expect_identical(rownames(X)[pairs[,2]], msubbed$gene2)
    
        msub.names <- paste0(msubbed$gene1, ".", msubbed$gene2)
        msub.ref <- ref.com[match(msub.names, ref.names),]
        msub.ref$FDR <- p.adjust(msub.ref$p.value, method="BH")
        rownames(msub.ref) <- NULL
        expect_equal(msub.ref, msubbed)
    
        pairs2 <- rownames(X)[pairs]
        dim(pairs2) <- dim(pairs)
        msubbed2 <- correlatePairs(X, pairings=pairs2)
        expect_equal(msubbed, msubbed2)
    
        # Checking for proper interaction with subset.row
        keep <- 1:10
        lsub2 <- correlatePairs(X, pairings=pairs, subset.row=keep)
        repairs <- pairs[pairs[,1] %in% keep & pairs[,2] %in% keep,,drop=FALSE]
        lexpected <- correlatePairs(X, pairings=repairs)
        expect_equal(lexpected, lsub2)
    }

    expect_error(correlatePairs(X, pairings=matrix(0, 0, 0)), 
        "'pairings' should be a numeric/character matrix with 2 columns")
    expect_error(correlatePairs(X, pairings=cbind(1,2,3)), 
        "'pairings' should be a numeric/character matrix with 2 columns")
    expect_error(correlatePairs(X, pairings=cbind(TRUE, FALSE)), 
        "'pairings' should be a numeric/character matrix with 2 columns")
})

####################################################################################################

set.seed(10003)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

test_that("correlatePairs works correctly with SingleCellExperiment objects", {
    ref <- correlatePairs(X)
    X2 <- SingleCellExperiment(list(logcounts=X))
    out <- correlatePairs(X2)
    expect_equal(out, ref)
})

test_that("correlatePairs fails properly upon silly inputs", {
    expect_identical(nrow(correlatePairs(X[0,])), 0L)
    expect_identical(nrow(correlatePairs(X[1,,drop=FALSE])), 0L)

    expect_true(all(is.nan(correlatePairs(X[,0])$rho)))
    expect_true(all(is.na(correlatePairs(X, design=diag(ncol(X)))$rho)))
    expect_true(all(is.nan(correlatePairs(X, block=1:ncol(X))$rho)))
})
