# Tests the pairwiseTTests function.
# library(scran); library(testthat); source("test-pairwise-t.R")

REFFUN <- function(y, grouping, direction="any", lfc=0) 
# A reference function using the t.test function.
{ 
    output <- pairwiseTTests(y, grouping, direction=direction, lfc=lfc)
    grouping <- factor(grouping)
    clust.vals <- levels(grouping)
    alt.hyp <- switch(direction, any="two.sided", up="greater", down="less")

    for (host in clust.vals) {
        host.y <- y[,grouping==host,drop=FALSE]
        for (target in setdiff(clust.vals, host)) {
            target.y <- y[,grouping==target,drop=FALSE]

            if (ncol(host.y)>1L && ncol(target.y)>1L) {
                pval <- numeric(nrow(y))
                for (i in seq_along(pval)) {

                    if (lfc==0) {
                        cur.p <- t.test(host.y[i,], target.y[i,], alternative=alt.hyp)$p.value
                    } else {
                        if (direction=="any") {
                            left.p <- t.test(host.y[i,], target.y[i,], alternative="less", mu=-lfc)$p.value + 
                                t.test(host.y[i,], target.y[i,], alternative="less", mu=lfc)$p.value
                            left.p <- left.p / 2
                            right.p <- t.test(host.y[i,], target.y[i,], alternative="greater", mu=-lfc)$p.value + 
                                t.test(host.y[i,], target.y[i,], alternative="greater", mu=lfc)$p.value
                            right.p <- right.p / 2
                            cur.p <- pmin(left.p, right.p) * 2
                        } else if (direction=="up") {
                            cur.p <- t.test(host.y[i,], target.y[i,], alternative=alt.hyp, mu=lfc)$p.value 
                        } else {
                            cur.p <- t.test(host.y[i,], target.y[i,], alternative=alt.hyp, mu=-lfc)$p.value 
                        }
                    }
                    pval[i] <- cur.p
                }
            } else {
                pval <- rep(NA_real_, nrow(y))
            }
            
			currow <- which(output$pairs[,1]==host & output$pairs[,2]==target)
            curres <- output$statistics[[currow]]
			expect_equal(unname(curres$logFC), unname(rowMeans(host.y) - rowMeans(target.y)))
            expect_equal(pval, curres$p.value)
            expect_equal(p.adjust(pval, method="BH"), curres$FDR)
            expect_identical(rownames(y), rownames(curres))
        }
    }  
    return(TRUE)
}

set.seed(70000)
ncells <- 200
ngenes <- 250
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
X <- scater::normalizeCounts(dummy, colSums(dummy))
rownames(X) <- seq_len(nrow(X))

set.seed(7000001)
test_that("pairwiseTTests works as expected without blocking or design matrices", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)

    REFFUN(X, clusters)
    REFFUN(X, clusters, direction="up")
    REFFUN(X, clusters, direction="down")

    REFFUN(X, clusters, lfc=0.2)
    REFFUN(X, clusters, lfc=0.2, direction="up")
    REFFUN(X, clusters, lfc=0.2, direction="down")

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

    out <- pairwiseTTests(X, re.clusters)
    ref <- pairwiseTTests(X, clusters)
    subset <- match(paste0(ref$pairs$first, ".", ref$pairs$second), 
        paste0(out$pairs$first, ".", out$pairs$second))
    expect_false(any(is.na(subset)))
    expect_equal(out$statistics[subset], ref$statistics)
})

###################################################################

BLOCKFUN <- function(y, grouping, block, direction="any", ...) {
    out <- pairwiseTTests(y, grouping, block=block, direction=direction, ...)
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
            subgroup <- grouping[chosen]

            N1 <- sum(subgroup==curpair[1])
            N2 <- sum(subgroup==curpair[2])
            if (N1==0 || N2==0) {
                next
            } 

            block.weights[[B]] <- 1/(1/N1 + 1/N2)
            block.res <- pairwiseTTests(y[,chosen], grouping[chosen], direction=direction, ...)
            to.use <- which(block.res$pairs$first==curpair[1] & block.res$pairs$second==curpair[2])

            block.lfc[[B]] <- block.res$statistics[[to.use]]$logFC
            if (direction=="any") { 
                # Recovering one-sided p-values for separate combining across blocks.
                going.up <- block.lfc[[B]] > 0
                curp <- block.res$statistics[[to.use]]$p.value/2
                block.up[[B]] <- ifelse(going.up, curp, 1-curp)
                block.down[[B]] <- ifelse(!going.up, curp, 1-curp)
            } else {
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

        # Doing the same for the p-values (for direction="any", we combine in each direction first).
        up.pval <- do.call(rbind, block.up)
        keep <- !is.na(up.pval[,1])
        if (any(keep)) { 
            up.Q <- qnorm(up.pval[keep,,drop=FALSE])
            up.Q <- colSums(up.Q * block.weights[keep]) / sqrt(sum(block.weights[keep]^2))
            up.p <- pnorm(up.Q)
            
            down.pval <- do.call(rbind, block.down)
            down.Q <- qnorm(down.pval[keep,,drop=FALSE])
            down.Q <- colSums(down.Q * block.weights[keep]) / sqrt(sum(block.weights[keep]^2))
            down.p <- pnorm(down.Q)
        } else {
            up.p <- down.p <- rep(NA_real_, nrow(ref.res))
        }

        if (direction=="any") {
            expect_equal(pmin(up.p, down.p) * 2, ref.res$p.value)
        } else if (direction=="up") {
            expect_equal(up.p, ref.res$p.value)
        } else if (direction=="down") {
            expect_equal(down.p, ref.res$p.value)
        }
    }

    return(TRUE)
}

set.seed(7000002)
test_that("pairwiseTTests works as expected with blocking", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)
    block <- sample(3, ncol(X), replace=TRUE)

    BLOCKFUN(X, clusters, block)
    BLOCKFUN(X, clusters, block, direction="up")
    BLOCKFUN(X, clusters, block, direction="down")

    BLOCKFUN(X, clusters, block, lfc=0.2)
    BLOCKFUN(X, clusters, block, lfc=0.2, direction="up")
    BLOCKFUN(X, clusters, block, lfc=0.2, direction="down")

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
    BLOCKFUN(X, re.clust, re.block)
})

###################################################################

LINEARFUN <- function(y, grouping, design, direction="any", lfc=0) {
    output <- pairwiseTTests(y, grouping, design=design, direction=direction, lfc=lfc)

    grouping <- factor(grouping)
    clust.vals <- levels(grouping)
    design2 <- model.matrix(~ 0 + grouping)
    colnames(design2) <- clust.vals
    design2 <- cbind(design2, design) # assume 'design' does not have an intercept.

    for (host in clust.vals) {
        design.custom <- design2
        design.custom[,host] <- 1
        fit <- limma::lmFit(y, design.custom)

        for (target in setdiff(clust.vals, host)) {
			currow <- which(output$pairs[,1]==host & output$pairs[,2]==target)
            curres <- output$statistics[[currow]]
            cur.lfc <- -fit$coefficients[,target] # Minus, as 'host' is currently the intercept.
			expect_equal(unname(curres$logFC), unname(cur.lfc))

            left <- pt((cur.lfc + lfc) / (fit$sigma * fit$stdev.unscaled[,target]), lower.tail=TRUE, df = fit$df.residual)
            right <- pt((cur.lfc - lfc) / (fit$sigma * fit$stdev.unscaled[,target]), lower.tail=FALSE, df = fit$df.residual)

            if (direction=="any") {
                left2 <- pt((cur.lfc - lfc) / (fit$sigma * fit$stdev.unscaled[,target]), lower.tail=TRUE, df = fit$df.residual)
                right2 <- pt((cur.lfc + lfc) / (fit$sigma * fit$stdev.unscaled[,target]), lower.tail=FALSE, df = fit$df.residual)
                left.p <- (left + left2)/2
                right.p <- (right + right2)/2
                pval <- pmin(left.p, right.p) * 2

            } else if (direction=="up") {
                pval <- right
            } else {
                pval <- left
            }

            pval <- unname(pval)
            expect_equal(pval, curres$p.value)
            expect_equal(p.adjust(pval, method="BH"), curres$FDR)
            expect_identical(rownames(y), rownames(curres))
        }
    }  
    return(TRUE)
}

set.seed(7000003)
test_that("pairwiseTTests works as expected with a design matrix", {
    clust <- kmeans(t(X), centers=3)
    clusters <- as.factor(clust$cluster)

    covariate <- cbind(runif(ncol(X)))
    LINEARFUN(X, clusters, covariate)
    LINEARFUN(X, clusters, covariate, direction="up")
    LINEARFUN(X, clusters, covariate, direction="down")

    alternative <- cbind(runif(ncol(X)), sample(0:1, ncol(X), replace=TRUE))
    LINEARFUN(X, clusters, alternative, lfc=0.2)
    LINEARFUN(X, clusters, alternative, lfc=0.2, direction="up")
    LINEARFUN(X, clusters, alternative, lfc=0.2, direction="down")
})

###################################################################

set.seed(7000004)
test_that("pairwiseTTests behaves as expected with subsetting", {
    y <- matrix(rnorm(12000), ncol=12)
    g <- gl(4,3)
    X <- cbind(runif(ncol(y)))

    # Integer subsetting.
    expect_identical(
        pairwiseTTests(y, g, subset.row=1:10),
        pairwiseTTests(y[1:10,], g)
    )
    expect_identical(
        pairwiseTTests(y, g, design=X, subset.row=1:10),
        pairwiseTTests(y[1:10,], g, design=X)
    )

    # Logical subsetting.
    keep <- rbinom(nrow(y), 1, 0.5)==1
    expect_identical(
        pairwiseTTests(y, g, subset.row=keep),
        pairwiseTTests(y[keep,], g, gene.names=which(keep))
    )
    expect_identical(
        pairwiseTTests(y, g, design=X, subset.row=keep),
        pairwiseTTests(y[keep,], g, design=X, gene.names=which(keep))
    )

    # Character subsetting.
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    chosen <- sample(rownames(y), 100)
    expect_identical(
        pairwiseTTests(y, g, subset.row=chosen),
        pairwiseTTests(y[chosen,], g)
    )
    expect_identical(
        pairwiseTTests(y, g, design=X, subset.row=chosen),
        pairwiseTTests(y[chosen,], g, design=X)
    )
})

set.seed(7000005)
test_that("pairwiseTTests behaves as expected with log-transformation", {
    y <- matrix(rnorm(12000), ncol=20)
    g <- gl(5,4)
    X <- cbind(rnorm(ncol(y)))

    # For Welch:
    ref <- pairwiseTTests(y, g)
    out <- pairwiseTTests(y, g, log.p=TRUE)
    expect_identical(ref$pairs, out$pairs)

    for (i in seq_along(ref$statistics)) {
        expect_equal(ref$statistics[[i]]$logFC, out$statistics[[i]]$logFC)
        expect_equal(log(ref$statistics[[i]]$p.value), out$statistics[[i]]$log.p.value)
        expect_equal(log(ref$statistics[[i]]$FDR), out$statistics[[i]]$log.FDR)
    }

    # For linear modelling:
    ref <- pairwiseTTests(y, g, design=X)
    out <- pairwiseTTests(y, g, design=X, log.p=TRUE)
    expect_identical(ref$pairs, out$pairs)

    for (i in seq_along(ref$statistics)) {
        expect_equal(ref$statistics[[i]]$logFC, out$statistics[[i]]$logFC)
        expect_equal(log(ref$statistics[[i]]$p.value), out$statistics[[i]]$log.p.value)
        expect_equal(log(ref$statistics[[i]]$FDR), out$statistics[[i]]$log.FDR)
    }
})

set.seed(7000006)
test_that("pairwiseTTests fails gracefully with silly inputs", {
    y <- matrix(rnorm(12000), ncol=20)
    g <- gl(5,4)
    X <- cbind(rnorm(ncol(y)))

    # Errors on incorrect inputs.
    expect_error(pairwiseTTests(y[,0], g), "does not equal")
    expect_error(pairwiseTTests(y, rep(1, ncol(y))), "need at least two")
    expect_error(pairwiseTTests(y, g, design=X[0,,drop=FALSE]), "is not equal")
    expect_error(pairwiseTTests(y, g, design=cbind(rep(1, ncol(y)))), "not of full rank")
    expect_error(pairwiseTTests(y, g, gene.names="A"), "not equal to the number of rows")

    # No genes.
    empty <- pairwiseTTests(y[0,], g)
    expect_identical(length(empty$statistics), nrow(empty$pairs))
    expect_true(all(sapply(empty$statistics, nrow)==0L))

    empty <- pairwiseTTests(y[0,], g, design=X)
    expect_identical(length(empty$statistics), nrow(empty$pairs))
    expect_true(all(sapply(empty$statistics, nrow)==0L))

    # Avoid NA p-values when variance is zero.
    clusters <- rep(1:2, each=ncol(y)/2)
    stuff <- matrix(clusters, ngenes, ncol(y), byrow=TRUE)
    out <- pairwiseTTests(stuff, clusters)
    expect_true(all(out$statistics[[1]]$FDR < 1e-8))
    expect_true(all(out$statistics[[2]]$FDR < 1e-8))
    expect_equal(out$statistics[[1]]$logFC, rep(-1, ngenes))
    expect_equal(out$statistics[[2]]$logFC, rep(1, ngenes))

    out <- pairwiseTTests(stuff, clusters, design=X)
    expect_true(all(out$statistics[[1]]$FDR < 1e-8))
    expect_true(all(out$statistics[[2]]$FDR < 1e-8))
    expect_equal(out$statistics[[1]]$logFC, rep(-1, ngenes))
    expect_equal(out$statistics[[2]]$logFC, rep(1, ngenes))
})
