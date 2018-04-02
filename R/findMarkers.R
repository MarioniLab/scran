#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
.findMarkers <- function(x, clusters, block=NULL, design=NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), 
                         lfc=0, log.p=FALSE, full.stats=FALSE, subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors in 'design' or 'block'.
#
# written by Aaron Lun
# created 22 March 2017
{
    ncells <- ncol(x)
    clusters <- as.factor(clusters)
    if (length(clusters)!=ncells) {
        stop("length of 'clusters' does not equal 'ncol(x)'")
    }
    pval.type <- match.arg(pval.type) 
    direction <- match.arg(direction)  
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Estimating the parameters.
    if (!is.null(block) || is.null(design)) {
        fit <- .test_block_internal(x, subset.row, clusters, block=block, direction=direction, lfc=lfc, full.stats=full.stats)
    } else {
        fit <- .fit_lm_internal(x, subset.row, clusters, design=design, direction=direction, lfc=lfc, full.stats=full.stats)
    }

    # Figuring out how to rank within each DataFrame.
    all.log.pval <- fit$p.value
    all.stats <- fit$stats
    output <- vector("list", length(all.stats))
    names(output) <- names(all.stats)

    for (host in names(output)) {
        stat.df <- do.call(DataFrame, c(all.stats[[host]], list(check.names=FALSE)))
        if (!full.stats) {
            colnames(stat.df) <- sprintf("logFC.%s", colnames(stat.df))
        } else {
            colnames(stat.df) <- sprintf("stats.%s", colnames(stat.df))
        }

        # Computing the combined p-value somehow.
        log.p.val <- do.call(cbind, all.log.pval[[host]])
        pval <- .combine_pvalues(log.p.val, pval.type=pval.type, log.p.in=TRUE, log.p.out=log.p)

        if (pval.type=="any") {
            # Ranking by position within each gene list.
            rank.out <- .rank_top_genes(log.p.val)
            min.rank <- rank.out$rank
            min.p <- rank.out$value
            gene.order <- order(min.rank, min.p)
            preamble <- DataFrame(Top=min.rank)
        } else {
            gene.order <- order(pval)
            if (log.p) {
                preamble <- DataFrame(log.IUT.p=pval)
            } else {
                preamble <- DataFrame(IUT.p=pval)
            }
        }

        # Performing a log-transformed version of the FDR correction, if desired.
        if (log.p) {
            preamble$log.FDR <- .logBH(pval)
        } else {
            preamble$FDR <- p.adjust(pval, method="BH")
        }

        # Running through all stats lists and computing the FDR.
        if (full.stats) {
            if (log.p) {
                for (target in names(stat.df)) {
                    stat.df[[target]]$log.FDR <- .logBH(stat.df[[target]]$log.p.value)
                }
            } else {
                for (target in names(stat.df)) {
                    cur.p <- exp(stat.df[[target]]$log.p.value)
                    stat.df[[target]]$p.value <- cur.p
                    stat.df[[target]]$log.p.value <- NULL
                    stat.df[[target]]$FDR <- p.adjust(cur.p, method="BH")
                }
            } 
        }

        # Producing the output object.
        gene.names <- rownames(x)[subset.row]
        if (is.null(gene.names)) { 
            gene.names <- subset.row
        }
        marker.set <- DataFrame(preamble, stat.df, check.names=FALSE, row.names=gene.names)
        marker.set <- marker.set[gene.order,,drop=FALSE]
        output[[host]] <- marker.set
    }

    return(output)
}

###########################################################
# Internal functions (blocking)
###########################################################

.test_block_internal <- function(x, subset.row, clusters, block=NULL, direction="any", lfc=0, full.stats=FALSE) 
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
        
    # Setting up block-level outputs. 
    out.means <- vector("list", length(by.block))
    out.s2 <- out.n <- out.means
    clust.vals <- levels(clusters)

    for (b in seq_along(by.block)) {
        curblock <- by.block[[b]]       
        by.cluster <- split(curblock - 1L, clusters[curblock], drop=FALSE)
        stopifnot(identical(names(by.cluster), clust.vals))
     
        # Calculating the statistics for each block.
        stats <- .Call(cxx_fit_oneway, by.cluster, x, subset.row - 1L)
        means <- stats[[1]]
        sigma2 <- stats[[2]]
        colnames(means) <- colnames(sigma2) <- clust.vals

        out.means[[b]] <- means
        out.s2[[b]] <- sigma2
        out.n[[b]] <- lengths(by.cluster)
    }

    # Setting up the output containers.
    out.p <- out.stats <- .create_output_container(clust.vals) 
    nblocks <- length(by.block)

    # Running through all pairs of comparisons.
    for (i in seq_along(clust.vals)) {
        host <- clust.vals[i]
        targets <- clust.vals[seq_len(i-1L)]

        for (target in targets) {
            all.lfc <- all.left <- all.right <- vector("list", nblocks)
            all.weight <- numeric(nblocks)
            valid.test <- logical(nblocks) 

            # Performing the same pairwise t-test across all blocks.
            for (b in seq_len(nblocks)) { 
                host.n <- out.n[[b]][host]
                target.n <- out.n[[b]][target]
                t.out <- .get_t_test_stats(host.s2=out.s2[[b]][,host], target.s2=out.s2[[b]][,target],
                                           host.n=host.n, target.n=target.n)
                
                cur.err <- t.out$err
                cur.df <- t.out$test.df
                cur.lfc <- out.means[[b]][,host] - out.means[[b]][,target]
                p.out <- .run_t_test(cur.lfc, cur.err, cur.df, thresh.lfc=lfc, direction=direction)

                all.lfc[[b]] <- cur.lfc
                all.left[[b]] <- p.out$left
                all.right[[b]] <- p.out$right

                # Weights are inversely proportional to the squared error of the log-fold change,
                # assuming equal variance across blocks and groups for simplicity.
                # (Using the actual variance is dangerous as some blocks have zero variance
                # with a defined log-fold change, if there are many zeroes.)
                all.weight[b] <- 1/(1/host.n + 1/target.n)

                # Indicating that there is no valid test statistic if the df is zero.
                valid.test[b] <- all(!is.na(cur.df))
            }

            # Combining the p-values and log-fold changes across blocks.
            com.left <- .run_stouffer(all.left, all.weight, valid.test)
            com.right <- .run_stouffer(all.right, all.weight, valid.test)

            # Flipping left/right to get the p-value from the reversed comparison.
            hvt.p <- .choose_leftright_pvalues(com.left, com.right, direction=direction) 
            out.p[[host]][[target]] <- hvt.p 
            tvh.p <- .choose_leftright_pvalues(com.right, com.left, direction=direction)
            out.p[[target]][[host]] <- tvh.p

            # Symmetrical log-fold changes, hence the -1.
            com.lfc <- .weighted_average_vals(all.lfc, all.weight, weighted=TRUE)
            if (!full.stats) { 
                out.stats[[host]][[target]] <- com.lfc
                out.stats[[target]][[host]] <- -com.lfc
            } else {
                # Alternatively, saving all of the stats as a DataFrame.
                out.stats[[host]][[target]] <- .create_full_stats(com.lfc, hvt.p, rownames(x))
                out.stats[[target]][[host]] <- .create_full_stats(-com.lfc, tvh.p, rownames(x))
            }
        }
    }
    return(list(p.value=out.p, stats=out.stats))    
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
    final <- .weighted_average_vals(all.Z, weights, weighted=TRUE)
    pnorm(final, log.p=TRUE)
}

###########################################################
# Internal functions (linear modelling)
###########################################################

#' @importFrom stats model.matrix
#' @importFrom limma lmFit contrasts.fit
.fit_lm_internal <- function(x, subset.row, clusters, design, direction="any", lfc=0, full.stats=FALSE)
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

    # Setting up output containers.
    out.p <- out.stats <- .create_output_container(clust.vals) 
    ngenes <- length(subset.row)

    # Doing a dummy fit, to avoid having to manually calculate standard errors.
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)
    for (h in seq_along(clust.vals)) { 
        host <- clust.vals[h]
        
        # Computing standard errors (via limma).
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
            out.p[[host]][[target]] <- hvt.p
            tvh.p <- .choose_leftright_pvalues(test.out$right, test.out$left, direction=direction)
            out.p[[target]][[host]] <- tvh.p 

            if (!full.stats) {
                out.stats[[host]][[target]] <- cur.lfc
                out.stats[[target]][[host]] <- -cur.lfc
            } else {
                out.stats[[host]][[target]] <- .create_full_stats(cur.lfc, hvt.p, rownames(x))
                out.stats[[target]][[host]] <- .create_full_stats(-cur.lfc, tvh.p, rownames(x))
            }
        }
    }
    return(list(p.value=out.p, stats=out.stats))
}

###########################################################
# Internal functions (p-value calculations)
###########################################################

#' @importFrom stats pt
.run_t_test <- function(cur.lfc, cur.err, cur.df, thresh.lfc=0, direction="any") 
# This runs the t-test given the relevant statistics. We use the TREAT method
# when thresh.lfc!=0, which tests against a null where the log-fold change is
# takes values at the extremes of [-thresh, thresh]. This is 50% distributed
# across both extremes, hence the -log(2) at the end.
{
    thresh.lfc <- abs(thresh.lfc)
    if (thresh.lfc==0) {    
        cur.t <- cur.lfc/sqrt(cur.err)
        left <- pt(cur.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
        right <- pt(cur.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
    } else {
        upper.t <- (cur.lfc - thresh.lfc)/cur.err
        lower.t <- (cur.lfc + thresh.lfc)/cur.err

        if (direction=="any") { 
            left <- .add_log_values(pt(upper.t, df=cur.df, lower.tail=TRUE, log.p=TRUE),
                                    pt(lower.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)) - log(2)
            right <- .add_log_values(pt(upper.t, df=cur.df, lower.tail=FALSE, log.p=TRUE),
                                     pt(lower.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)) - log(2)
        } else {
            # Note that if direction='up', only 'right' is used; nonetheless, 'left' is still calculated
            # using 'lower.t' so as to allow quick calculation of the p-value for the reversed contrast,
            # by simply swapping 'left' and 'right'. The same applies when direction='down'.
            left <- pt(lower.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
            right <- pt(upper.t, df=cur.df, lower.tail=FALSE, log.p=TRUE) 
        }
    }
    return(list(left=left, right=right))
}

.add_log_values <- function(x, y) 
# This performs log(exp(x) + exp(y)) in a numerically
# stable fashion, i.e., without actually running the 'exp's.
{
    pmax(x, y) + log(1 + exp(-abs(x-y)))
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

.combine_pvalues <- function(P, pval.type="any", log.p.in=FALSE, log.p.out=FALSE) 
# This function combines the p-values using Simes' method or via the IUT.
# Additional arguments are involved to specify whether the input/output should be logged.
{
    if (pval.type=="any") { 
        # Computing the Simes p-value (with NA protection).
        P <- .Call(cxx_combine_simes, P, log.p.in)
    } else {
        # Computing the IUT p-value.
        P <- P[.find_largest_col(P)]
    }

    # Deciding what to return (at this point, P is the same log-status as it was supplied).
    if (log.p.in && !log.p.out) { 
        P <- exp(P) 
    } else if (!log.p.in && log.p.out) {
        P <- log(P)
    }
    return(P)
}

.logBH <- function(log.p.val) 
# Same as log(p.adjust(exp(log.p.val), method="BH")), without
# the need to undo and redo the log-transformations.
{
    o <- order(log.p.val)
    repval <- log.p.val[o] + log(length(o)/seq_along(o))
    repval <- rev(cummin(rev(repval)))
    repval[o] <- repval
    return(repval)
}

###########################################################
# Internal functions (other)
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

.rank_top_genes <- function(metrics) 
# This computes the rank and the minimum metric for each gene.
{
    ngenes <- nrow(metrics)
    ncon <- ncol(metrics)
    min.rank <- min.val <- rep(NA_integer_, ngenes)

    for (con in seq_len(ncon)) { 
        cur.val <- metrics[,con]
        cur.rank <- rank(cur.val, ties.method="first", na.last="keep")
        min.rank <- pmin(min.rank, cur.rank, na.rm=TRUE)
        min.val <- pmin(min.val, cur.val, na.rm=TRUE)
    }
    
    return(list(rank=min.rank, value=min.val))
}

.find_largest_col <- function(metrics) {
    metrics[is.na(metrics)] <- -Inf # protection against NAs.
    largest <- max.col(metrics)
    cbind(seq_along(largest), largest)
}

#' @importFrom S4Vectors DataFrame
.create_full_stats <- function(cur.lfc, cur.p, gene.names) {
    I(DataFrame(logFC=cur.lfc, log.p.value=as.vector(cur.p), 
                check.names=FALSE, row.names=gene.names))
}

###########################################################
# S4 method definitions
###########################################################

#' @export
setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

#' @export
setMethod("findMarkers", "ANY", .findMarkers)

#' @importFrom SummarizedExperiment assay
#' @export
setMethod("findMarkers", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .findMarkers(assay(x, i=assay.type), ..., subset.row=subset.row)
})                                 


