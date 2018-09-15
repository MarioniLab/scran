#' @export
#' @importFrom S4Vectors DataFrame
pairwiseTTests <- function(x, clusters, block=NULL, design=NULL, direction=c("any", "up", "down"),
    lfc=0, log.p=FALSE, gene.names=rownames(x), subset.row=NULL)
# Performs pairwise Welch t-tests between clusters.
#
# written by Aaron Lun
# created 15 September 2018
{
    ncells <- ncol(x)
    clusters <- as.factor(clusters)
    if (length(clusters)!=ncells) {
        stop("length of 'clusters' does not equal 'ncol(x)'")
    }
    if (nlevels(clusters) < 2L) {
        stop("need at least two unique levels in 'clusters'")
    }

    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    if (is.null(gene.names)) {
        gene.names <- subset.row
    } else if (length(gene.names)!=nrow(x)) {
        stop("length of 'gene.names' is not equal to the number of rows")
    } else {
        gene.names <- gene.names[subset.row]
    }

    direction <- match.arg(direction)

    if (!is.null(block) || is.null(design)) {
        results <- .test_block_internal(x, subset.row, clusters, block=block, direction=direction, lfc=lfc, gene.names=gene.names, log.p=log.p)
    } else {
        results <- .fit_lm_internal(x, subset.row, clusters, design=design, direction=direction, lfc=lfc, gene.names=gene.names, log.p=log.p)
    }

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names))
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

###########################################################
# Internal functions (blocking)
###########################################################

#' @importFrom S4Vectors DataFrame
.test_block_internal <- function(x, subset.row, clusters, block=NULL, direction="any", lfc=0, gene.names=NULL, log.p=TRUE)
# This looks at every level of the blocking factor and performs
# t-tests between pairs of clusters within each blocking level.
{
    ncells <- ncol(x)
    if (is.null(block)) {
        by.block <- list(`1`=seq_len(ncells))
    } else {
        if (length(block)!=ncells) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        by.block <- split(seq_len(ncells), block)
    }

    # Calculating the statistics for each block.
    out.means <- vector("list", length(by.block))
    out.s2 <- out.n <- out.means
    clust.vals <- levels(clusters)

    for (b in seq_along(by.block)) {
        curblock <- by.block[[b]]
        by.cluster <- split(curblock - 1L, clusters[curblock], drop=FALSE)
        stopifnot(identical(names(by.cluster), clust.vals))

        stats <- .Call(cxx_fit_oneway, by.cluster, x, subset.row - 1L)
        means <- stats[[1]]
        sigma2 <- stats[[2]]
        colnames(means) <- colnames(sigma2) <- clust.vals

        out.means[[b]] <- means
        out.s2[[b]] <- sigma2
        out.n[[b]] <- lengths(by.cluster)
    }

    # Running through all pairs of comparisons.
    out.stats <- .create_output_container(clust.vals)
    nblocks <- length(by.block)

    for (i in seq_along(clust.vals)) {
        host <- clust.vals[i]
        targets <- clust.vals[seq_len(i-1L)]

        for (target in targets) {
            all.lfc <- all.left <- all.right <- vector("list", nblocks)
            all.weight <- numeric(nblocks)
            valid.test <- logical(nblocks)

            # Performing the same pairwise t-test within each block.
            for (b in seq_len(nblocks)) {
                host.n <- out.n[[b]][host]
                target.n <- out.n[[b]][target]
                t.out <- .get_t_test_stats(host.s2=out.s2[[b]][,host], target.s2=out.s2[[b]][,target], host.n=host.n, target.n=target.n)

                cur.err <- t.out$err
                cur.df <- t.out$test.df
                cur.lfc <- out.means[[b]][,host] - out.means[[b]][,target]
                p.out <- .run_t_test(cur.lfc, cur.err, cur.df, thresh.lfc=lfc, direction=direction)

                all.lfc[[b]] <- cur.lfc
                all.left[[b]] <- p.out$left
                all.right[[b]] <- p.out$right

                # Weights are inversely proportional to the squared error of the log-fold change,
                # _assuming equal variance_ across blocks and groups for simplicity.
                # (Using the actual variance is dangerous as some blocks have zero variance
                # with a defined log-fold change, if there are many zeroes.)
                all.weight[b] <- 1/(1/host.n + 1/target.n)

                # Indicating that there is no valid test statistic if the df is zero.
                valid.test[b] <- all(!is.na(cur.df))
            }

            # Combining the p-values for each side across blocks.
            com.left <- .run_stouffer(all.left, all.weight, valid.test)
            com.right <- .run_stouffer(all.right, all.weight, valid.test)

            # Flipping left/right to get the p-value from the reversed comparison.
            hvt.p <- .choose_leftright_pvalues(com.left, com.right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(com.right, com.left, direction=direction)

            # Symmetrical log-fold changes, hence the -1.
            com.lfc <- .weighted_average_vals(all.lfc, all.weight, weighted=TRUE)
            out.stats[[host]][[target]] <- .create_full_stats(com.lfc, hvt.p, gene.names, log.p=log.p)
            out.stats[[target]][[host]] <- .create_full_stats(-com.lfc, tvh.p, gene.names, log.p=log.p)
        }
    }
	
    out.stats
}

.get_t_test_stats <- function(host.s2, target.s2, host.n, target.n)
# Computes the error variance of the log-fold change and the
# degrees of freedom for the t-distribution.
{
    host.df <- max(0L, host.n - 1L)
    target.df <- max(0L, target.n - 1L)

    # Avoid unlikely but potential problems with discreteness.
    host.s2 <- pmax(host.s2, 1e-8)
    target.s2 <- pmax(target.s2, 1e-8)

    if (host.df > 0L && target.df > 0L) {
        # Perform Welch's t-test here.
        host.err <- host.s2/host.n
        target.err <- target.s2/target.n
        cur.err <- host.err + target.err
        cur.df <- cur.err^2 / (host.err^2/host.df + target.err^2/target.df)
    } else {
        cur.err <- cur.df <- NA_real_
    }
    return(list(err=cur.err, test.df=cur.df))
}

#' @importFrom stats qnorm pnorm
.run_stouffer <- function(log.pvals, weights, valid)
# Uses Stouffer's method to combine the one-sided p-values, with weights for each column.
# Only valid tests are used, i.e., with non-zero d.f. for computing the test statistic.
# We also skip the calculations if there's only one or no tests.
{
    ngenes <- length(log.pvals[[1]])

    if (!all(valid)) {
        log.pvals <- log.pvals[valid]
        weights <- weights[valid]
    }

    NC <- length(log.pvals)
    if (NC==1L) {
        return(log.pvals[[1]])
    } else if (NC==0L) {
        return(rep(NA_real_, ngenes))
    }

    all.Z <- lapply(log.pvals, qnorm, log.p=TRUE)
    final <- Reduce("+", mapply(all.Z, weights, FUN="*", SIMPLIFY=FALSE)) / sqrt(sum(weights^2))
    pnorm(final, log.p=TRUE)
}

###########################################################
# Internal functions (linear modelling)
###########################################################

#' @importFrom stats model.matrix
#' @importFrom limma lmFit contrasts.fit
.fit_lm_internal <- function(x, subset.row, clusters, design, direction="any", lfc=0, gene.names=NULL, log.p=TRUE)
# This fits a linear model to each gene and performs a t-test for
# differential expression between clusters.
{
    full.design <- model.matrix(~0 + clusters)
    clust.vals <- levels(clusters)
    colnames(full.design) <- clust.vals

    # Removing terms to avoid linear dependencies on the intercept.
    if (nrow(design)!=ncol(x)) {
        stop("'nrow(design)' is not equal to 'ncol(x)'")
    }
    out <- qr.solve(design, cbind(rep(1, nrow(design))))
    to.drop <- abs(out) > 1e-8
    if (any(to.drop)) {
        design <- design[,-which(to.drop)[1],drop=FALSE]
    }
    full.design <- cbind(full.design, design) # Other linear dependencies will trigger errors in .ranksafe_QR.

    # Getting coefficient estimates (need to undo column pivoting to get the actual estimates).
    QR <- .ranksafe_qr(full.design)
    resid.df <- nrow(full.design) - ncol(full.design)
    if (resid.df <= 0L) {
        stop("no residual d.f. in design matrix for variance estimation")
    }
    stats <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row - 1L, TRUE)
    coefficients <- stats[[1]]
    coefficients[QR$pivot,] <- coefficients
    sigma2 <- stats[[3]]

    # Running through every pair of clusters.
    out.stats <- .create_output_container(clust.vals)
    ngenes <- length(subset.row)
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)

    for (h in seq_along(clust.vals)) {
        host <- clust.vals[h]

        # Computing standard errors (via limma, to avoid having to manually calculate standard errors).
        con <- matrix(0, ncol(full.design), h-1L)
        diag(con) <- -1
        con[h,] <- 1
        lfit2 <- contrasts.fit(lfit, con)

        # Computing log-fold changes, t-statistics and p-values on a per-contrast basis.
        # This _could_ be vectorised, but it's less confusing to do it like this,
        # and there's not much speed gain to be had from vectorizing over contrasts.
        ref.coef <- coefficients[h,]
        for (tdex in seq_len(h-1L)) {
            target <- clust.vals[tdex]
            cur.lfc <- ref.coef - coefficients[tdex,]

            test.out <- .run_t_test(cur.lfc, lfit2$stdev.unscaled[tdex]*sqrt(sigma2), resid.df, thresh.lfc=lfc, direction=direction)
            hvt.p <- .choose_leftright_pvalues(test.out$left, test.out$right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(test.out$right, test.out$left, direction=direction)

            out.stats[[host]][[target]] <- .create_full_stats(cur.lfc, hvt.p, gene.names, log.p=log.p)
            out.stats[[target]][[host]] <- .create_full_stats(-cur.lfc, tvh.p, gene.names, log.p=log.p)
        }
    }

    out.stats
}

###########################################################
# Internal functions (t-test calculations)
###########################################################

#' @importFrom stats pt
.run_t_test <- function(cur.lfc, cur.err, cur.df, thresh.lfc=0, direction="any")
# This runs the t-test given the relevant statistics, regardless of how
# they were computed (i.e., within blocks, or in a linear model).
{
    thresh.lfc <- abs(thresh.lfc)
    if (thresh.lfc==0) {
        cur.t <- cur.lfc/sqrt(cur.err)
        left <- pt(cur.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
        right <- pt(cur.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
    } else {
        upper.t <- (cur.lfc - thresh.lfc)/sqrt(cur.err)
        lower.t <- (cur.lfc + thresh.lfc)/sqrt(cur.err)

        left.lower <- pt(lower.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
        right.upper <- pt(upper.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)

        if (direction=="any") {
            # Using the TREAT method, which tests against a null where the log-fold change is
            # takes values at the extremes of [-thresh, thresh]. The null probability is 50%
            # distributed across both extremes, hence the -log(2) at the end.
            left.upper <- pt(upper.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
            right.lower <- pt(lower.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
            left <- .add_log_values(left.upper, left.lower) - log(2)
            right <- .add_log_values(right.upper, right.lower) - log(2)
        } else {
            # For one-sided tests, testing against the lfc threshold in the specified direction.
            # Note that if direction='up', only 'right' is used; nonetheless, 'left' is still calculated
            # using 'lower.t' so as to allow quick calculation of the p-value for the reversed contrast,
            # by simply swapping 'left' and 'right'. The same applies when direction='down'.
            left <- left.lower
            right <- right.upper
        }
    }
    return(list(left=left, right=right))
}

.add_log_values <- function(x, y)
# This performs log(exp(x) + exp(y)) in a numerically
# stable fashion, i.e., without actually running the 'exp's.
{
    pmax(x, y) + log1p(exp(-abs(x-y)))
}

.choose_leftright_pvalues <- function(left, right, direction="any")
# Choosing whether to use the p-values from the left (lower.tail), or
# the p-values from the right (upper.tail), or to make it two-sided.
# Assumes 'left' and 'right' are log-transformed p-values.
{
    if (direction=="up") {
        return(right)
    } else if (direction=="down") {
        return(left)
    } else {
        log.p.out <- pmin(left, right) + log(2)
        return(log.p.out)
    }
}

###########################################################
# Internal functions (output formatting)
###########################################################

.create_output_container <- function(clust.vals) {
    out <- vector("list", length(clust.vals))
    names(out) <- clust.vals
    for (i in seq_along(clust.vals)) {
        targets <- clust.vals[-i]
        collected <- vector("list", length(targets))
        names(collected) <- targets
        host <- clust.vals[i]
        out[[host]] <- collected
    }
    return(out)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
.create_full_stats <- function(cur.lfc, cur.p, gene.names, log.p=TRUE) {
    cur.p <- as.vector(cur.p)
    if (log.p) {
        DataFrame(logFC=cur.lfc, log.p.value=cur.p, log.FDR=.logBH(cur.p), check.names=FALSE, row.names=gene.names)
    } else {
        cur.p <- exp(cur.p)
        DataFrame(logFC=cur.lfc, p.value=cur.p, FDR=p.adjust(cur.p, method="BH"), check.names=FALSE, row.names=gene.names)
    }
}
