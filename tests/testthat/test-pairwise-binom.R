# Tests the pairwiseBinom() function.
# library(scran); library(testthat); source("setup.R"); source("test-pairwise-binom.R")

REFFUN <- function(y, grouping, direction="any", lfc=0) 
# A reference function using the t.test function.
{ 
    output <- pairwiseBinom(y, grouping, direction=direction, lfc=lfc)
    grouping <- factor(grouping)
    clust.vals <- levels(grouping)
    alt.hyp <- switch(direction, any="two.sided", up="greater", down="less")

    for (host in clust.vals) {
        host.i <- grouping==host
        host.y <- rowSums(y[,host.i,drop=FALSE] > 0)
        host.n <- sum(host.i)

        for (target in setdiff(clust.vals, host)) {
            target.i <- grouping==target
            target.y <- rowSums(y[,target.i,drop=FALSE] > 0)
            target.n <- sum(target.i)

            n <- host.y + target.y
            p <- host.n/(host.n + target.n)

            if (host.n > 0L && target.n > 0L) {
                lvals <- edgeR::cpm(cbind(host.y, target.y), lib.size=c(host.n, target.n), 
                    prior.count=1, log=TRUE)
                effect <- unname(lvals[,1] - lvals[,2])

                pval <- numeric(nrow(y))
                for (i in seq_along(pval)) {
                    cur.n <- n[i]
                    cur.y <- host.y[i]

                    if (cur.n) {
                        if (direction=="any") {
                            right <- binom.test(cur.y, cur.n, p=p, alternative="greater")
                            left <- binom.test(cur.y, cur.n, p=p, alternative="less")
                            pv <- pmin(right$p.value, left$p.value, 0.5) * 2
                        } else {
                            pv <- binom.test(cur.y, cur.n, p=p, alternative=alt.hyp)$p.value
                        }
                    } else {
                        pv <- 1
                    }
                    pval[i] <- unname(pv)
                }
            } else {
                pval <- effect <- rep(NA_real_, nrow(y))
            }

			currow <- which(output$pairs[,1]==host & output$pairs[,2]==target)
            curres <- output$statistics[[currow]]
			expect_equal(unname(curres$logFC), effect)
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
means <- runif(ngenes, 0, 5)
X <- matrix(rpois(ngenes*ncells, lambda=means), ncol=ncells, nrow=ngenes)
rownames(X) <- seq_len(nrow(X))

set.seed(8000001)
test_that("pairwiseBinom works as expected without blocking", {
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
    REFFUN(X, re.clust)

    # Checking what happens if there is an empty level.
    re.clusters <- clusters
    levels(re.clusters) <- 1:4

    expect_warning(out <- pairwiseBinom(X, re.clusters), "no within-block")
    ref <- pairwiseBinom(X, clusters)
    subset <- match(paste0(ref$pairs$first, ".", ref$pairs$second), 
        paste0(out$pairs$first, ".", out$pairs$second))
    expect_false(any(is.na(subset)))
    expect_equal(out$statistics[subset], ref$statistics)
})

set.seed(800000112)
test_that("pairwiseBinom works as expected with a log-fold change threshold", {
    # Throwing in a very small lfc to check that 
    # the fundamental calculations are executed correctly
    # for the lfc-based function.
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)

    out <- pairwiseBinom(X, clusters)
    ref <- pairwiseBinom(X, clusters, lfc=1e-8)
    expect_equal(out, ref, tol=1e-6)

    out <- pairwiseBinom(X, clusters, direction="up")
    ref <- pairwiseBinom(X, clusters, lfc=1e-8, direction="up")
    expect_equal(out, ref, tol=1e-6)

    out <- pairwiseBinom(X, clusters, direction="down")
    ref <- pairwiseBinom(X, clusters, lfc=1e-8, direction="down")
    expect_equal(out, ref, tol=1e-6)

    # Just getting some test coverage here, not much that can be done
    # without rewriting all of the relevant code for 'p'.
    out <- pairwiseBinom(X, clusters, lfc=0.5)
    out <- pairwiseBinom(X, clusters, lfc=0.5, direction="up")
    out <- pairwiseBinom(X, clusters, lfc=0.5, direction="down")
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
test_that("pairwiseBinom responds to non-standard level ordering", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)
    f1 <- factor(clusters)
    f2 <- factor(clusters, rev(levels(f1)))
    FACTORCHECK(pairwiseBinom(X, f1), pairwiseBinom(X, f2))
})

set.seed(80000012)
test_that("pairwiseBinom responds to restriction", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)

    restrict <- c("B", "C")
    keep <- clusters %in% restrict
    expect_identical(pairwiseBinom(X, clusters, restrict=restrict),
       pairwiseBinom(X[,keep], clusters[keep]))

    restrict <- c("A", "D", "E")
    keep <- clusters %in% restrict
    expect_identical(pairwiseBinom(X, clusters, restrict=restrict),
       pairwiseBinom(X[,keep], clusters[keep]))

    exclude <- c("A", "B", "C")
    keep <- !clusters %in% exclude
    expect_identical(pairwiseBinom(X, clusters, exclude=exclude),
       pairwiseBinom(X[,keep], clusters[keep]))
})

###################################################################

BLOCKFUN <- function(y, grouping, block, direction="any", ...) {
    out <- pairwiseBinom(y, grouping, block=block, direction=direction, ...)
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
            subgroup <- factor(grouping[chosen]) # refactoring to get rid of empty levels.

            N1 <- sum(subgroup==curpair[1])
            N2 <- sum(subgroup==curpair[2])
            if (N1==0 || N2==0) {
                next
            } 
            block.weights[[B]] <- N1 + N2

            suby <- y[,chosen,drop=FALSE]
            if (direction=="any") { 
                # Recovering one-sided p-values for separate combining across blocks.
                block.res.up <- pairwiseBinom(suby, subgroup, direction="up", ...)
                to.use.up <- which(block.res.up$pairs$first==curpair[1] & block.res.up$pairs$second==curpair[2])
                block.res.down <- pairwiseBinom(suby, subgroup, direction="down", ...)
                to.use.down <- which(block.res.down$pairs$first==curpair[1] & block.res.down$pairs$second==curpair[2])

                block.lfc[[B]] <- block.res.up$statistics[[to.use.up]]$logFC
                block.up[[B]] <- block.res.up$statistics[[to.use.up]]$p.value
                block.down[[B]] <- block.res.down$statistics[[to.use.down]]$p.value
            } else {
                block.res <- pairwiseBinom(suby, subgroup, direction=direction, ...)
                to.use <- which(block.res$pairs$first==curpair[1] & block.res$pairs$second==curpair[2])
                block.lfc[[B]] <- block.res$statistics[[to.use]]$logFC
                block.up[[B]] <- block.down[[B]] <- block.res$statistics[[to.use]]$p.value
            }
        }

        block.weights <- unlist(block.weights)
        if (length(block.weights)==0) {
            expect_equal(ref.res$logFC, rep(NA_real_, nrow(ref.res)))
            expect_equal(ref.res$p.value, rep(NA_real_, nrow(ref.res)))
            next
        }

        # Taking a weighted average.
        all.lfc <- do.call(rbind, block.lfc)
        ave.lfc <- colSums(all.lfc * block.weights) / sum(block.weights)
        expect_equal(ave.lfc, ref.res$logFC)

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
test_that("pairwiseBinom works as expected with blocking", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)
    block <- sample(3, ncol(X), replace=TRUE)

    BLOCKFUN(X, clusters, block)
    BLOCKFUN(X, clusters, block, direction="up")
    BLOCKFUN(X, clusters, block, direction="down")

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
test_that("pairwiseBinom with blocking works across multiple cores", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)
    block <- sample(3, ncol(X), replace=TRUE)
    ref <- pairwiseBinom(X, clusters, block=block)

    expect_equal(ref, pairwiseBinom(X, clusters, block=block, BPPARAM=safeBPParam(2)))

    expect_equal(ref, pairwiseBinom(X, clusters, block=block, BPPARAM=SnowParam(2)))
})

set.seed(80000022)
test_that("pairwiseBinom with blocking responds to non-standard level ordering", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)
    f1 <- factor(clusters)
    f2 <- factor(clusters, rev(levels(f1)))

    b <- sample(1:3, ncol(X), replace=TRUE)
    FACTORCHECK(pairwiseBinom(X, f1, block=b), pairwiseBinom(X, f2, block=b))

    b1 <- factor(b, 1:3)
    b2 <- factor(b, 3:1)
    FACTORCHECK(pairwiseBinom(X, f1, block=b1), pairwiseBinom(X, f2, block=b2))
})

set.seed(80000023)
test_that("pairwiseBinom with blocking responds to restriction", {
    clusters <- sample(LETTERS[1:5], ncol(X), replace=TRUE)

    restrict <- c("B", "C")
    keep <- clusters %in% restrict
    b <- sample(1:3, ncol(X), replace=TRUE)
    expect_identical(pairwiseBinom(X, clusters, restrict=restrict, block=b),
       pairwiseBinom(X[,keep], clusters[keep], block=b[keep]))

    restrict <- c("A", "D", "E")
    keep <- clusters %in% restrict
    expect_identical(pairwiseBinom(X, clusters, restrict=restrict, block=b),
       pairwiseBinom(X[,keep], clusters[keep], block=b[keep]))

    # What happens if the block and cluster are correlated?
    b2 <- b
    b2[!clusters %in% restrict] <- 0
    expect_identical(pairwiseBinom(X, clusters, restrict=restrict, block=b2),
       pairwiseBinom(X[,keep], clusters[keep], block=b2[keep]))
})

###################################################################

set.seed(8000004)
test_that("pairwiseBinom behaves as expected with subsetting", {
    y <- matrix(rnorm(1200), ncol=12)
    rownames(y) <- seq_len(nrow(y))
    g <- gl(4,3)
    X <- cbind(runif(ncol(y)))

    # Integer subsetting.
    expect_identical(
        pairwiseBinom(y, g, subset.row=1:10),
        pairwiseBinom(y[1:10,], g)
    )

    # Logical subsetting.
    keep <- rbinom(nrow(y), 1, 0.5)==1
    expect_identical(
        pairwiseBinom(y, g, subset.row=keep),
        pairwiseBinom(y[keep,], g) 
    )

    # Character subsetting.
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    chosen <- sample(rownames(y), 100)
    expect_identical(
        pairwiseBinom(y, g, subset.row=chosen),
        pairwiseBinom(y[chosen,], g)
    )

    # Auto-generates names for the subset.
    y <- y0 <- matrix(rnorm(1200), ncol=12)
    rownames(y) <- seq_len(nrow(y))
    expect_identical(
        pairwiseBinom(y0, g, subset.row=10:1),
        pairwiseBinom(y[10:1,], g)
    )
})

set.seed(8000005)
test_that("pairwiseBinom behaves as expected with log-transformation", {
    y <- matrix(rnorm(1200), ncol=20)
    g <- gl(5,4)
    X <- cbind(rnorm(ncol(y)))

    # For Welch:
    ref <- pairwiseBinom(y, g)
    out <- pairwiseBinom(y, g, log.p=TRUE)
    expect_identical(ref$pairs, out$pairs)

    for (i in seq_along(ref$statistics)) {
        expect_equal(ref$statistics[[i]]$effect, out$statistics[[i]]$effect)
        expect_equal(log(ref$statistics[[i]]$p.value), out$statistics[[i]]$log.p.value)
        expect_equal(log(ref$statistics[[i]]$FDR), out$statistics[[i]]$log.FDR)
    }
})

set.seed(80000051)
test_that("pairwiseBinom works with SEs and SCEs", {
    y <- matrix(rnorm(1200), ncol=12)
    g <- gl(4,3)

    out <- pairwiseBinom(y, g)
    out2 <- pairwiseBinom(SummarizedExperiment(list(logcounts=y)), g)
    expect_identical(out, out2)

    X2 <- SingleCellExperiment(list(logcounts=y))
    colLabels(X2) <- g
    out3 <- pairwiseBinom(X2)
    expect_identical(out, out3)
})

set.seed(8000006)
test_that("pairwiseBinom fails gracefully with silly inputs", {
    y <- matrix(rnorm(1200), ncol=20)
    g <- gl(5,4)

    # Errors on incorrect inputs.
    expect_error(pairwiseBinom(y[,0], g), "does not equal")
    expect_error(pairwiseBinom(y, rep(1, ncol(y))), "need at least two")

    # No genes.
    empty <- pairwiseBinom(y[0,], g)
    expect_identical(length(empty$statistics), nrow(empty$pairs))
    expect_true(all(sapply(empty$statistics, nrow)==0L))

    # Avoid NA p-values when variance is zero.
    stuff <- matrix(c(0, 1), ngenes, ncol(y))
    out <- pairwiseBinom(stuff, g)
    expect_true(all(out$statistics[[1]]$FDR==1))
    expect_true(all(out$statistics[[2]]$FDR==1))
    expect_equal(out$statistics[[1]]$logFC, rep(0, ngenes))
    expect_equal(out$statistics[[2]]$logFC, rep(0, ngenes))
})
