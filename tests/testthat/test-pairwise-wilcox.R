# Tests the pairwiseWilcox() function.
# library(scran); library(testthat); source("setup.R"); source("test-pairwise-wilcox.R")

REFFUN <- function(y, grouping, direction="any", lfc=0) 
# A reference function using the t.test function.
{ 
    output <- pairwiseWilcox(y, grouping, direction=direction, lfc=lfc)
    grouping <- factor(grouping)
    clust.vals <- levels(grouping)
    alt.hyp <- switch(direction, any="two.sided", up="greater", down="less")

    for (host in clust.vals) {
        host.y <- y[,grouping==host,drop=FALSE]
        for (target in setdiff(clust.vals, host)) {
            target.y <- y[,grouping==target,drop=FALSE]

            if (ncol(host.y) * ncol(target.y) > 1L) {
                auc <- pval <- numeric(nrow(y))
                if (lfc==0) { 
                    for (i in seq_along(pval)) {
                        result <- wilcox.test(host.y[i,], target.y[i,], alternative=alt.hyp, exact=FALSE)
                        auc[i] <- result$statistic 
                        pval[i] <- result$p.value
                    }
                } else {
                    for (i in seq_along(pval)) {
                        host.vals <- host.y[i,]
                        target.vals <- target.y[i,]
                        if (direction=="any") {
                            left.result1 <- wilcox.test(host.vals, target.vals, alternative="less", mu=-lfc, exact=FALSE)
                            left.result2 <- wilcox.test(host.vals, target.vals, alternative="less", mu=lfc, exact=FALSE)
                            right.p <- wilcox.test(host.vals, target.vals, alternative="greater", mu=lfc, exact=FALSE)$p.value
                            auc[i] <- (left.result1$statistic + left.result2$statistic) / 2
                            pval[i] <- pmin(left.result1$p.value, right.p, 0.5) * 2
                        } else if (direction=="up") {
                            result <- wilcox.test(host.vals, target.vals, alternative=alt.hyp, mu=lfc, exact=FALSE)
                            auc[i] <- result$statistic
                            pval[i] <- result$p.value
                        } else {
                            result <- wilcox.test(host.vals, target.vals, alternative=alt.hyp, mu=-lfc, exact=FALSE)
                            auc[i] <- result$statistic 
                            pval[i] <- result$p.value
                        }
                    }
                }
                auc <- auc / (ncol(host.y) * ncol(target.y))
            } else {
                pval <- auc <- rep(NA_real_, nrow(y))
            }
            
			currow <- which(output$pairs[,1]==host & output$pairs[,2]==target)
            curres <- output$statistics[[currow]]
			expect_equal(unname(curres$AUC), auc)
            expect_equal(pval, curres$p.value)
            expect_equal(p.adjust(pval, method="BH"), curres$FDR)
            expect_identical(rownames(y), rownames(curres))
        }
    }  
    return(TRUE)
}

set.seed(80000)
ncells <- 200
ngenes <- 150
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
X <- scuttle::normalizeCounts(dummy, colSums(dummy))
rownames(X) <- seq_len(nrow(X))

set.seed(8000001)
test_that("pairwiseWilcox works as expected without blocking", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)

    REFFUN(X, clusters)
    REFFUN(X, clusters, direction="up")
    REFFUN(X, clusters, direction="down")

    # Checking what happens if one of the groups has only one element.
    re.clust <- clust$cluster
    re.clust[1] <- 4
    re.clust <- factor(re.clust)
    REFFUN(X, re.clust)

    # Checking what happens if two of the groups have only one element.
    re.clust <- clust$cluster
    re.clust[1:2] <- 4:5
    re.clust <- factor(re.clust)
    expect_warning(REFFUN(X, re.clust), "no within-block")

    # Checking what happens if there is an empty level.
    re.clusters <- clusters
    levels(re.clusters) <- 1:4

    expect_warning(out <- pairwiseWilcox(X, re.clusters), "no within-block")
    ref <- pairwiseWilcox(X, clusters)
    subset <- match(paste0(ref$pairs$first, ".", ref$pairs$second), 
        paste0(out$pairs$first, ".", out$pairs$second))
    expect_false(any(is.na(subset)))
    expect_equal(out$statistics[subset], ref$statistics)
})

set.seed(80000011)
test_that("pairwiseWilcox works with a log-fold change threshold", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)

    # Note that some numerical imprecision means that the p-values 
    # may be slightly off, depending on the clustering.
    REFFUN(X, clusters, lfc=0.5)
    REFFUN(X, clusters, direction="up", lfc=0.5)
    REFFUN(X, clusters, direction="down", lfc=0.5)
})

FACTORCHECK <- function(left, right) {
    expect_identical(names(left), names(right))

    oL <- order(left$pairs[,1], left$pairs[,2])
    oR <- order(right$pairs[,1], right$pairs[,2])
    expect_identical(left$pairs[oL,], right$pairs[oR,])

    expect_identical(names(left$statistics)[oL], names(right$statistics)[oR])
    for (x in seq_along(oL)) {
        curleft <- left$statistics[[oL[x]]]
        curright <- right$statistics[[oR[x]]]
        expect_identical(sort(colnames(curleft)), sort(colnames(curright)))
        expect_equal(curleft, curright[,colnames(curleft)])
    }
    return(TRUE)
}

set.seed(80000011)
test_that("pairwiseWilcox responds to non-standard level ordering", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)
    f1 <- factor(clusters)
    f2 <- factor(clusters, rev(levels(f1)))
    FACTORCHECK(pairwiseWilcox(X, f1), pairwiseWilcox(X, f2))
})

set.seed(80000012)
test_that("pairwiseWilcox responds to restriction", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)

    restrict <- c("B", "C")
    keep <- clusters %in% restrict
    expect_identical(pairwiseWilcox(X, clusters, restrict=restrict),
       pairwiseWilcox(X[,keep], clusters[keep]))

    restrict <- c("A", "D", "E")
    keep <- clusters %in% restrict
    expect_identical(pairwiseWilcox(X, clusters, restrict=restrict),
       pairwiseWilcox(X[,keep], clusters[keep]))

    exclude <- c("A", "B", "C")
    keep <- !clusters %in% exclude
    expect_identical(pairwiseWilcox(X, clusters, exclude=exclude),
       pairwiseWilcox(X[,keep], clusters[keep]))
})

set.seed(80000013)
test_that("pairwiseWilcox works correctly with lots of ties and zeros", {
    dummy <- matrix(rpois(ngenes*ncells, lambda=1), ncol=ncells, nrow=ngenes)
    rownames(dummy) <- seq_len(nrow(dummy))
    clusters <- factor(sample(4, ncells, replace=TRUE))

    REFFUN(dummy, clusters)
    REFFUN(dummy, clusters, direction="up")
    REFFUN(dummy, clusters, direction="down")
    REFFUN(dummy, clusters, lfc=0.5)

    # Same behavior with the sparse matrix.
    sparse <- as(dummy, "dgCMatrix")
    expect_identical(pairwiseWilcox(dummy, clusters), pairwiseWilcox(sparse, clusters))
    expect_identical(pairwiseWilcox(dummy, clusters, lfc=1.7), pairwiseWilcox(sparse, clusters, lfc=1.7))
})

###################################################################

BLOCKFUN <- function(y, grouping, block, direction="any", ...) {
    out <- pairwiseWilcox(y, grouping, block=block, direction=direction, ...)
    ngroups <- length(unique(grouping))
    expect_equal(nrow(out$pairs), ngroups^2L - ngroups)
    expect_identical(nrow(out$pairs), length(out$statistics))

    for (p in seq_len(nrow(out$pairs))) {
        curpair <- unlist(out$pairs[p,])
        ref.res <- out$statistics[[p]]

        # Extracting block-wise results.
        block.weights <- block.up <- block.down <- block.lfc <- list()
        for (b in unique(block)) { 
            B <- as.character(b)
            chosen <- block==b & grouping %in% curpair
            subgroup <- factor(grouping[chosen]) # refactoring to eliminate unused levels.

            N1 <- sum(subgroup==curpair[1])
            N2 <- sum(subgroup==curpair[2])
            if (N1==0 || N2==0) {
                next
            } 
            block.weights[[B]] <- N1 * N2

            suby <- y[,chosen,drop=FALSE]
            if (direction=="any") { 
                # Recovering one-sided p-values for separate combining across blocks.
                block.res.up <- pairwiseWilcox(suby, subgroup, direction="up", ...)
                to.use.up <- which(block.res.up$pairs$first==curpair[1] & block.res.up$pairs$second==curpair[2])
                block.up[[B]] <- block.res.up$statistics[[to.use.up]]$p.value

                block.res.down <- pairwiseWilcox(suby, subgroup, direction="down", ...)
                to.use.down <- which(block.res.down$pairs$first==curpair[1] & block.res.down$pairs$second==curpair[2])
                block.down[[B]] <- block.res.down$statistics[[to.use.down]]$p.value

                block.lfc[[B]] <- block.res.up$statistics[[to.use.up]]$AUC
            } else {
                block.res <- pairwiseWilcox(suby, subgroup, direction=direction, ...)
                to.use <- which(block.res$pairs$first==curpair[1] & block.res$pairs$second==curpair[2])
                block.lfc[[B]] <- block.res$statistics[[to.use]]$AUC
                block.up[[B]] <- block.down[[B]] <- block.res$statistics[[to.use]]$p.value
            }
        }

        block.weights <- unlist(block.weights)
        if (length(block.weights)==0) {
            expect_equal(ref.res$AUC, rep(NA_real_, nrow(ref.res)))
            expect_equal(ref.res$p.value, rep(NA_real_, nrow(ref.res)))
            next
        }

        # Taking a weighted average.
        all.lfc <- do.call(rbind, block.lfc)
        ave.lfc <- colSums(all.lfc * block.weights) / sum(block.weights)
        expect_equal(ave.lfc, ref.res$AUC)

        # Combining p-values in each direction.
        up.p <- metapod::parallelStouffer(block.up, weights=block.weights)$p.value
        down.p <- metapod::parallelStouffer(block.down, weights=block.weights)$p.value

        if (direction=="any") {
            expect_equal(pmin(up.p, down.p, 0.5) * 2, ref.res$p.value)
        } else if (direction=="up") {
            expect_equal(up.p, ref.res$p.value)
        } else if (direction=="down") {
            expect_equal(down.p, ref.res$p.value)
        }
    }

    return(TRUE)
}

set.seed(8000002)
test_that("pairwiseWilcox works as expected with blocking", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.character(clust$cluster)
    block <- sample(3, ncol(X), replace=TRUE)

    BLOCKFUN(X, clusters, block)
    BLOCKFUN(X, clusters, block, direction="up")
    BLOCKFUN(X, clusters, block, direction="down")

    # Not checking direction='any', because we can't easily recover one-sided
    # p-values when we need to account for the discreteness of the test statistic.
    BLOCKFUN(X, clusters, block, direction="up", lfc=0.5)
    BLOCKFUN(X, clusters, block, direction="down", lfc=0.5)

    # Checking what happens to a block-specific group.
    re.clust <- clust$cluster
    re.clust[block!=1 & re.clust==1] <- 2
    re.clust <- factor(re.clust)
    BLOCKFUN(X, re.clust, block)

    # Checking what happens to a group-specific block.
    re.clust <- clust$cluster
    re.clust[block==1] <- 1
    re.clust <- factor(re.clust)
    BLOCKFUN(X, re.clust, block)

    # Checking what happens to a doubly-specific group and block.
    re.clust <- clust$cluster
    re.clust[block==1] <- 1
    re.block <- block
    re.block[re.clust==1] <- 1
    expect_warning(BLOCKFUN(X, re.clust, re.block), "no within-block")
})

set.seed(80000021)
test_that("pairwiseWilcox with blocking works across multiple cores", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)
    block <- sample(3, ncol(X), replace=TRUE)
    ref <- pairwiseWilcox(X, clusters, block=block)

    expect_equal(ref, pairwiseWilcox(X, clusters, block=block, BPPARAM=safeBPParam(2)))

    expect_equal(ref, pairwiseWilcox(X, clusters, block=block, BPPARAM=SnowParam(2)))
})

set.seed(80000022)
test_that("pairwiseWilcox with blocking responds to non-standard level ordering", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)
    f1 <- factor(clusters)
    f2 <- factor(clusters, rev(levels(f1)))

    b <- sample(1:3, ncol(X), replace=TRUE)
    FACTORCHECK(pairwiseWilcox(X, f1, block=b), pairwiseWilcox(X, f2, block=b))

    b1 <- factor(b, 1:3)
    b2 <- factor(b, 3:1)
    FACTORCHECK(pairwiseWilcox(X, f1, block=b1), pairwiseWilcox(X, f2, block=b2))
})

set.seed(80000023)
test_that("pairwiseWilcox with blocking responds to restriction", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)

    restrict <- c("B", "C")
    keep <- clusters %in% restrict
    b <- sample(1:3, ncol(X), replace=TRUE)
    expect_identical(pairwiseWilcox(X, clusters, restrict=restrict, block=b),
       pairwiseWilcox(X[,keep], clusters[keep], block=b[keep]))

    restrict <- c("A", "D", "E")
    keep <- clusters %in% restrict
    expect_identical(pairwiseWilcox(X, clusters, restrict=restrict, block=b),
       pairwiseWilcox(X[,keep], clusters[keep], block=b[keep]))

    # What happens if the block and cluster are correlated?
    b2 <- b
    b2[!clusters %in% restrict] <- 0
    expect_identical(pairwiseWilcox(X, clusters, restrict=restrict, block=b2),
       pairwiseWilcox(X[,keep], clusters[keep], block=b2[keep]))
})

###################################################################

set.seed(8000004)
test_that("pairwiseWilcox behaves as expected with subsetting", {
    y <- matrix(rnorm(1200), ncol=12)
    rownames(y) <- seq_len(nrow(y))
    g <- gl(4,3)
    X <- cbind(runif(ncol(y)))

    # Integer subsetting.
    expect_identical(
        pairwiseWilcox(y, g, subset.row=1:10),
        pairwiseWilcox(y[1:10,], g)
    )

    # Logical subsetting.
    keep <- rbinom(nrow(y), 1, 0.5)==1
    expect_identical(
        pairwiseWilcox(y, g, subset.row=keep),
        pairwiseWilcox(y[keep,], g) 
    )

    # Character subsetting.
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    chosen <- sample(rownames(y), 100)
    expect_identical(
        pairwiseWilcox(y, g, subset.row=chosen),
        pairwiseWilcox(y[chosen,], g)
    )

    # Auto-generates names for the subset.
    y <- y0 <- matrix(rnorm(1200), ncol=12)
    rownames(y) <- seq_len(nrow(y))
    expect_identical(
        pairwiseWilcox(y0, g, subset.row=10:1),
        pairwiseWilcox(y[10:1,], g)
    )
})

set.seed(8000005)
test_that("pairwiseWilcox behaves as expected with log-transformation", {
    y <- matrix(rnorm(1200), ncol=20)
    g <- gl(5,4)
    X <- cbind(rnorm(ncol(y)))

    # For Welch:
    ref <- pairwiseWilcox(y, g)
    out <- pairwiseWilcox(y, g, log.p=TRUE)
    expect_identical(ref$pairs, out$pairs)

    for (i in seq_along(ref$statistics)) {
        expect_equal(ref$statistics[[i]]$AUC, out$statistics[[i]]$AUC)
        expect_equal(log(ref$statistics[[i]]$p.value), out$statistics[[i]]$log.p.value)
        expect_equal(log(ref$statistics[[i]]$FDR), out$statistics[[i]]$log.FDR)
    }
})

set.seed(80000051)
test_that("pairwiseWilcox works with SEs and SCEs", {
    y <- matrix(rnorm(1200), ncol=12)
    g <- gl(4,3)

    out <- pairwiseWilcox(y, g)
    out2 <- pairwiseWilcox(SummarizedExperiment(list(logcounts=y)), g)
    expect_identical(out, out2)

    X2 <- SingleCellExperiment(list(logcounts=y))
    colLabels(X2) <- g
    out3 <- pairwiseWilcox(X2)
    expect_identical(out, out3)
})

set.seed(8000006)
test_that("pairwiseWilcox fails gracefully with silly inputs", {
    y <- matrix(rnorm(1200), ncol=20)
    g <- gl(5,4)

    # Errors on incorrect inputs.
    expect_error(pairwiseWilcox(y[,0], g), "does not equal")
    expect_error(pairwiseWilcox(y, rep(1, ncol(y))), "need at least two")

    # No genes.
    empty <- pairwiseWilcox(y[0,], g)
    expect_identical(length(empty$statistics), nrow(empty$pairs))
    expect_true(all(sapply(empty$statistics, nrow)==0L))

    # Avoid NA p-values when variance is zero.
    clusters <- rep(1:2, each=ncol(y)/2)
    stuff <- matrix(clusters, ngenes, ncol(y), byrow=TRUE)
    out <- pairwiseWilcox(stuff, clusters)
    expect_true(all(out$statistics[[1]]$FDR < 1e-4))
    expect_true(all(out$statistics[[2]]$FDR < 1e-4))
    expect_equal(out$statistics[[1]]$AUC, rep(0, ngenes))
    expect_equal(out$statistics[[2]]$AUC, rep(1, ngenes))

    stuff <- matrix(0, ngenes, ncol(y))
    out <- pairwiseWilcox(stuff, clusters)
    expect_true(all(out$statistics[[1]]$p.value==1))
    expect_true(all(out$statistics[[2]]$p.value==1))
    expect_equal(out$statistics[[1]]$AUC, rep(0.5, ngenes))
    expect_equal(out$statistics[[2]]$AUC, rep(0.5, ngenes))
})
