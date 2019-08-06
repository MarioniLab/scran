#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
pairwiseTTests <- function(x, clusters, block=NULL, design=NULL, direction=c("any", "up", "down"),
    lfc=0, std.lfc=FALSE, log.p=FALSE, gene.names=rownames(x), subset.row=NULL, BPPARAM=SerialParam())
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
        results <- .test_block_internal(x, subset.row, clusters, block=block, direction=direction, lfc=lfc, 
            std.lfc=std.lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
    } else {
        results <- .fit_lm_internal(x, subset.row, clusters, design=design, direction=direction, lfc=lfc, 
            std.lfc=std.lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
    }

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names), use.names=FALSE)
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

###########################################################
# Internal functions (blocking)
###########################################################

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply SerialParam
.test_block_internal <- function(x, subset.row, clusters, block=NULL, direction="any", lfc=0, std.lfc=FALSE,
    gene.names=NULL, log.p=TRUE, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# t-tests between pairs of clusters within each blocking level.
{
    nclusters <- nlevels(clusters)
    if (is.null(block)) {
        all.clusters <- factor(as.integer(clusters), seq_len(nclusters))
        nblocks <- 1L
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- factor(block)
        nblocks <- nlevels(block)
        all.clusters <- factor(as.integer(clusters) + (as.integer(block) - 1L) * nclusters, seq_len(nblocks*nclusters))
    }

    # Calculating the statistics for each block.
    all.blocks <- split(seq_along(all.clusters) - 1L, all.clusters)
    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row, wout)
    by.core <- .split_matrix_by_workers(x, by.core)

    raw.stats <- bplapply(by.core, FUN=compute_blocked_stats_none, bygroup=all.blocks, BPPARAM=BPPARAM)
    all.means <- do.call(rbind, lapply(raw.stats, FUN=function(x) t(x[[1]])))
    all.vars <- do.call(rbind, lapply(raw.stats, FUN=function(x) t(x[[2]])))
    all.n <- table(all.clusters)

    clust.vals <- levels(clusters)
    out.s2 <- out.means <- out.n <- vector("list", nblocks)
    for (b in seq_len(nblocks)) {
        chosen <- (b-1L) * nclusters + seq_len(nclusters)
        means <- all.means[,chosen,drop=FALSE]
        sigma2 <- all.vars[,chosen,drop=FALSE]
        cur.n <- as.vector(all.n[chosen])

        colnames(means) <- colnames(sigma2) <- names(cur.n) <- clust.vals
        out.means[[b]] <- means
        out.s2[[b]] <- sigma2
        out.n[[b]] <- cur.n 
    }

    # Running through all pairs of comparisons.
    STATFUN <- function(b, host, target) {
        host.n <- out.n[[b]][host]
        target.n <- out.n[[b]][target]
        host.s2 <- out.s2[[b]][,host]
        target.s2 <- out.s2[[b]][,target]
        t.out <- .get_t_test_stats(host.s2=host.s2, target.s2=target.s2, host.n=host.n, target.n=target.n)

        cur.err <- t.out$err
        cur.df <- t.out$test.df
        cur.lfc <- out.means[[b]][,host] - out.means[[b]][,target]
        p.out <- .run_t_test(cur.lfc, cur.err, cur.df, thresh.lfc=lfc, direction=direction)

        effect.size <- cur.lfc
        if (std.lfc) {
            # Computing Cohen's D.
            pooled.s2 <- ((host.n - 1) * host.s2 + (target.n - 1) * target.s2)/(target.n + host.n - 2)
            effect.size <- effect.size / sqrt(pooled.s2)
        }
        
        list(effect=effect.size, left=p.out$left, right=p.out$right, 

            # Weights are inversely proportional to the squared error of the log-fold change,
            # _assuming equal variance_ across blocks and groups for simplicity.
            # (Using the actual variance is dangerous as some blocks have zero variance
            # with a defined log-fold change, if there are many zeroes.)
            weight=1/(1/host.n + 1/target.n),

            # Indicating that there is no valid test statistic if the df is zero.
            valid=all(!is.na(cur.df))
        )
    }

    .pairwise_blocked_template(x, clust.vals, nblocks, direction=direction,
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, FLIPFUN=function(x) -x, effect.name="logFC")
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

###########################################################
# Internal functions (linear modelling)
###########################################################

#' @importFrom stats model.matrix
#' @importFrom limma lmFit contrasts.fit
#' @importFrom BiocParallel bplapply SerialParam
.fit_lm_internal <- function(x, subset.row, clusters, design, direction="any", lfc=0, std.lfc=FALSE,
    gene.names=NULL, log.p=TRUE, BPPARAM=SerialParam())
# This fits a linear model to each gene and performs a t-test for
# differential expression between clusters.
{
    full.design <- model.matrix(~0 + clusters)
    clust.vals <- levels(clusters)
    colnames(full.design) <- clust.vals

    if (nrow(design)!=ncol(x)) {
        stop("'nrow(design)' is not equal to 'ncol(x)'")
    }
    full.design <- cbind(full.design, design) # Other linear dependencies will trigger errors in .ranksafe_QR.

    # Getting coefficient estimates (need to undo column pivoting to get the actual estimates).
    QR <- .ranksafe_qr(full.design)
    resid.df <- nrow(full.design) - ncol(full.design)
    if (resid.df <= 0L) {
        stop("no residual d.f. in design matrix for variance estimation")
    }

    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row-1L, wout)
    raw.stats <- bplapply(by.core, FUN=fit_linear_model, qr=QR$qr, qraux=QR$qraux, exprs=x, get_coefs=TRUE, BPPARAM=BPPARAM)
    coefficients <- do.call(cbind, lapply(raw.stats, "[[", i=1))
    coefficients[QR$pivot,] <- coefficients
    sigma2 <- unlist(lapply(raw.stats, "[[", i=3))
    sigma2 <- pmax(sigma2, 1e-8) # avoid unlikely but possible problems with discreteness.

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

            test.out <- .run_t_test(cur.lfc, lfit2$stdev.unscaled[tdex]^2*sigma2, resid.df, thresh.lfc=lfc, direction=direction)
            hvt.p <- .choose_leftright_pvalues(test.out$left, test.out$right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(test.out$right, test.out$left, direction=direction)

            if (std.lfc) {
                # Computing Cohen's D.
                cur.lfc <- cur.lfc / sqrt(sigma2)
            }
            out.stats[[host]][[target]] <- .create_full_stats(cur.lfc, p=hvt.p, gene.names=gene.names, log.p=log.p)
            out.stats[[target]][[host]] <- .create_full_stats(-cur.lfc, p=tvh.p, gene.names=gene.names, log.p=log.p)
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
