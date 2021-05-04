#' Perform pairwise t-tests
#' 
#' Perform pairwise Welch t-tests between groups of cells, possibly after blocking on uninteresting factors of variation.
#' 
#' @param x A numeric matrix-like object of normalized log-expression values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param groups A vector of length equal to \code{ncol(x)}, specifying the group assignment for each cell.
#' If \code{x} is a SingleCellExperiment, this is automatically derived from \code{\link{colLabels}}.
#' @param block A factor specifying the blocking level for each cell.
#' @param design A numeric matrix containing blocking terms for uninteresting factors.
#' Note that these factors should not be confounded with \code{groups}.
#' @param restrict A vector specifying the levels of \code{groups} for which to perform pairwise comparisons.
#' @param exclude A vector specifying the levels of \code{groups} for which \emph{not} to perform pairwise comparisons.
#' @param direction A string specifying the direction of log-fold changes to be considered in the alternative hypothesis.
#' @param lfc A positive numeric scalar specifying the log-fold change threshold to be tested against.
#' @param std.lfc A logical scalar indicating whether log-fold changes should be standardized.
#' @param log.p A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param gene.names Deprecated, use \code{row.data} in \code{\link{findMarkers}} instead to add custom annotation.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param assay.type A string specifying which assay values to use, usually \code{"logcounts"}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' 
#' @details
#' This function performs t-tests to identify differentially expressed genes (DEGs) between pairs of groups of cells.
#' A typical aim is to use the DEGs to determine cluster identity based on expression of marker genes with known biological activity.
#' A list of tables is returned where each table contains per-gene statistics for a comparison between one pair of groups.
#' This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.
#' 
#' We use t-tests as they are simple, fast and perform reasonably well for single-cell data (Soneson and Robinson, 2018).
#' However, if one of the groups contains fewer than two cells, no p-value will be reported for comparisons involving that group.
#' A warning will also be raised about insufficient degrees of freedom (d.f.) in such cases.
#' 
#' When \code{log.p=TRUE}, the log-transformed p-values and FDRs are reported using the natural base.
#' This is useful in cases with many cells such that reporting the p-values directly would lead to double-precision underflow.
#'
#' If \code{restrict} is specified, comparisons are only performed between pairs of groups in \code{restrict}.
#' This can be used to focus on DEGs distinguishing between a subset of the groups (e.g., closely related cell subtypes).
#'
#' If \code{exclude} is specified, comparisons are not performed between groups in \code{exclude}.
#' Similarly, if any entries of \code{groups} are \code{NA}, the corresponding cells are are ignored.
#' 
#' @section Direction and magnitude of the log-fold change:
#' Log-fold changes are reported as differences in the values of \code{x}.
#' Thus, all log-fold changes have the same base as whatever was used to perform the log-transformation in \code{x}.
#' If \code{\link{logNormCounts}} was used, this would be base 2.
#' 
#' If \code{direction="any"}, two-sided tests will be performed for each pairwise comparisons between groups.
#' Otherwise, one-sided tests in the specified direction will be used instead.
#' This can be used to focus on genes that are upregulated in each group of interest, which is often easier to interpret when assigning cell type to a cluster.
#' 
#' To interpret the setting of \code{direction}, consider the DataFrame for group X, in which we are comparing to another group Y.
#' If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in group X compared to Y.
#' If \code{direction="down"}, genes will only be significant if they are downregulated in group X compared to Y.
#' 
#' The magnitude of the log-fold changes can also be tested by setting \code{lfc}.
#' By default, \code{lfc=0} meaning that we will reject the null upon detecting any differential expression.
#' If this is set to some other positive value, the null hypothesis will change depending on \code{direction}:
#' \itemize{
#' \item If \code{direction="any"}, the null hypothesis is that the true log-fold change lies inside \code{[-lfc, lfc]}.
#' To be conservative, we perform one-sided tests against the boundaries of this interval, and combine them to obtain a two-sided p-value.
#' \item If \code{direction="up"}, the null hypothesis is that the true log-fold change is less than or equal to \code{lfc}.
#' A one-sided p-value is computed against the boundary of this interval.
#' \item If \code{direction="down"}, the null hypothesis is that the true log-fold change is greater than or equal to \code{-lfc}.
#' A one-sided p-value is computed against the boundary of this interval.
#' }
#' This is similar to the approach used in \code{\link[limma:eBayes]{treat}} and allows users to focus on genes with strong log-fold changes.
#' 
#' If \code{std.lfc=TRUE}, the log-fold change for each gene is standardized by the variance.
#' When the Welch t-test is being used, this is equivalent to Cohen's d.
#' Standardized log-fold changes may be more appealing for visualization as it avoids large fold changes due to large variance.
#' The choice of \code{std.lfc} does not affect the calculation of the p-values.
#' 
#' @section Blocking on uninteresting factors:
#' If \code{block} is specified, Welch t-tests are performed between groups within each level of \code{block}.
#' For each pair of groups, the p-values for each gene across all levels of \code{block} are combined using Stouffer's weighted Z-score method.
#' The reported log-fold change for each gene is also a weighted average of log-fold changes across levels.
#' 
#' The weight for a particular level is defined as \eqn{(1/N_x + 1/N_y)^{-1}}, 
#' where \eqn{Nx} and \eqn{Ny} are the number of cells in groups X and Y, respectively, for that level. 
#' This is inversely proportional to the expected variance of the log-fold change, provided that all groups and blocking levels have the same variance.
#' 
#' % In theory, a better weighting scheme would be to use the estimated standard error of the log-fold change to compute the weight.
#' % This would be more responsive to differences in variance between blocking levels, focusing on levels with low variance and high power.
#' % However, this is not safe in practice as genes with many zeroes can have very low standard errors, dominating the results inappropriately.
#' 
#' When comparing two groups, blocking levels are ignored if no p-value was reported, e.g., if there were insufficient cells for a group in a particular level. 
#' This includes levels that contain fewer than two cells for either group, as this cannot yield a p-value from the Welch t-test.
#' If all levels are ignored in this manner, the entire comparison will only contain \code{NA} p-values and a warning will be emitted.
#' 
#' @section Regressing out unwanted factors:
#' If \code{design} is specified, a linear model is instead fitted to the expression profile for each gene.
#' This linear model will include the \code{groups} as well as any blocking factors in \code{design}.
#' A t-test is then performed to identify DEGs between pairs of groups, using the values of the relevant coefficients and the gene-wise residual variance.
#'
#' Note that \code{design} must be full rank when combined with the \code{groups} terms, 
#' i.e., there should not be any confounding variables.
#' We make an exception for the common situation where \code{design} contains an \code{"(Intercept)"} column,
#' which is automatically detected and removed (emitting a warning along the way).
#' 
#' We recommend using \code{block} instead of \code{design} for uninteresting categorical factors of variation.
#' The former accommodates differences in the variance of expression in each group via Welch's t-test.
#' As a result, it is more robust to misspecification of the groups, as misspecified groups (and inflated variances) do not affect the inferences for other groups.
#' Use of \code{block} also avoids assuming additivity of effects between the blocking factors and the groups.
#' 
#' Nonetheless, use of \code{design} is unavoidable when blocking on real-valued covariates.
#' It is also useful for ensuring that log-fold changes/p-values are computed for comparisons between all pairs of groups
#' (assuming that \code{design} is not confounded with the groups).
#' This may not be the case with \code{block} if a pair of groups never co-occur in a single blocking level. 
#' 
#' @return
#' A list is returned containing \code{statistics} and \code{pairs}.
#' 
#' The \code{statistics} element is itself a list of \linkS4class{DataFrame}s.
#' Each DataFrame contains the statistics for a comparison between a pair of groups,
#' including the log-fold changes, p-values and false discovery rates.
#' 
#' The \code{pairs} element is a DataFrame where each row corresponds to an entry of \code{statistics}.
#' This contains the fields \code{first} and \code{second}, 
#' specifying the two groups under comparison in the corresponding DataFrame in \code{statistics}.
#' 
#' In each DataFrame in \code{statistics}, the log-fold change represents the change in the \code{first} group compared to the \code{second} group.
#' 
#' @author
#' Aaron Lun
#' 
#' @references
#' Whitlock MC (2005). 
#' Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
#' \emph{J. Evol. Biol.} 18, 5:1368-73.
#' 
#' Soneson C and Robinson MD (2018). 
#' Bias, robustness and scalability in single-cell differential expression analysis. 
#' \emph{Nat. Methods}
#' 
#' Lun ATL (2018).
#' Comments on marker detection in \emph{scran}.
#' \url{https://ltla.github.io/SingleCellThoughts/software/marker_detection/comments.html}
#'
#' @seealso
#' \code{\link{t.test}}, on which this function is based.
#'
#' \code{\link{combineMarkers}}, to combine pairwise comparisons into a single DataFrame per group.
#'
#' \code{\link{getTopMarkers}}, to obtain the top markers from each pairwise comparison.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(sce)), centers=3) 
#' 
#' # Vanilla application:
#' out <- pairwiseTTests(logcounts(sce), groups=kout$cluster)
#' out
#' 
#' # Directional with log-fold change threshold:
#' out <- pairwiseTTests(logcounts(sce), groups=kout$cluster, 
#'     direction="up", lfc=0.2)
#' out
#' 
#' @name pairwiseTTests
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle .subset2index
.pairwiseTTests <- function(x, groups, block=NULL, design=NULL, restrict=NULL, exclude=NULL, direction=c("any", "up", "down"),
    lfc=0, std.lfc=FALSE, log.p=FALSE, gene.names=NULL, subset.row=NULL, BPPARAM=SerialParam())
{
    groups <- .setup_groups(groups, x, restrict=restrict, exclude=exclude)
    direction <- match.arg(direction)

    if (!is.null(block) && !is.null(design)) {
        stop("cannot specify both 'block' and 'design'")
    } else if (!is.null(design)) {
        .fit_lm_internal(x, subset.row, groups, design=design, direction=direction, lfc=lfc, 
            std.lfc=std.lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
    } else {
        .test_block_internal(x, subset.row, groups, block=block, direction=direction, lfc=lfc, 
            std.lfc=std.lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
    }
}

#' @export
#' @rdname pairwiseTTests
setGeneric("pairwiseTTests", function(x, ...) standardGeneric("pairwiseTTests"))

#' @export
#' @rdname pairwiseTTests
setMethod("pairwiseTTests", "ANY", .pairwiseTTests)

#' @export
#' @rdname pairwiseTTests
#' @importFrom SummarizedExperiment assay
setMethod("pairwiseTTests", "SummarizedExperiment", function(x, ..., assay.type="logcounts") { 
    .pairwiseTTests(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname pairwiseTTests
#' @importFrom SingleCellExperiment colLabels
setMethod("pairwiseTTests", "SingleCellExperiment", function(x, groups=colLabels(x, onAbsence="error"), ...) {
    callNextMethod(x=x, groups=groups, ...)
})

###########################################################
# Internal functions (blocking)
###########################################################

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp 
.test_block_internal <- function(x, subset.row, groups, block=NULL, direction="any", lfc=0, std.lfc=FALSE,
    gene.names=NULL, log.p=TRUE, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# t-tests between pairs of groups within each blocking level.
{
    ngroups <- nlevels(groups)
    if (is.null(block)) {
        all.groups <- factor(as.integer(groups), seq_len(ngroups))
        nblocks <- 1L
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block[is.na(groups)] <- NA # avoid overstating 'nblocks'.
        block <- factor(block)
        nblocks <- nlevels(block)
        all.groups <- factor(as.integer(groups) + (as.integer(block) - 1L) * ngroups, seq_len(nblocks*ngroups))
    }

    # Calculating the statistics for each block.
    stats <- .compute_mean_var(x, BPPARAM=BPPARAM, subset.row=subset.row, design=NULL, 
        block.FUN=compute_blocked_stats_none, block=all.groups)

    all.means <- stats$means
    all.vars <- stats$vars
    all.n <- stats$ncells

    clust.vals <- levels(groups)
    out.s2 <- out.means <- out.n <- vector("list", nblocks)
    for (b in seq_len(nblocks)) {
        chosen <- (b-1L) * ngroups + seq_len(ngroups)
        means <- all.means[,chosen,drop=FALSE]
        sigma2 <- all.vars[,chosen,drop=FALSE]
        cur.n <- as.vector(all.n[chosen])

        dimnames(means) <- dimnames(sigma2) <- list(NULL, clust.vals)
        names(cur.n) <- clust.vals

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
        p.out <- .run_t_test(cur.lfc, cur.err, cur.df, thresh.lfc=lfc) 

        effect.size <- cur.lfc
        if (std.lfc) {
            # Computing Cohen's D.
            pooled.s2 <- ((host.n - 1) * host.s2 + (target.n - 1) * target.s2)/(target.n + host.n - 2)
            is.zero <- effect.size == 0
            effect.size <- effect.size / sqrt(pooled.s2)
            effect.size[is.zero] <- 0
        }
        
        list(forward=effect.size, reverse=-effect.size, left=p.out$left, right=p.out$right, 

            # Weights are inversely proportional to the squared error of the log-fold change,
            # _assuming equal variance_ across blocks and groups for simplicity.
            # (Using the actual variance is dangerous as some blocks have zero variance
            # with a defined log-fold change, if there are many zeroes.)
            weight=1/(1/host.n + 1/target.n),

            # Indicating that there is no valid test statistic if the df is zero.
            valid=all(!is.na(cur.df))
        )
    }

    gene.names <- .setup_gene_names(gene.names, x, subset.row)

    .pairwise_blocked_template(clust.vals, nblocks=nblocks, direction=direction,
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, effect.name="logFC",
        BPPARAM=BPPARAM)
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
#' @importFrom scuttle fitLinearModel
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame
.fit_lm_internal <- function(x, subset.row, groups, design, direction="any", lfc=0, std.lfc=FALSE,
    gene.names=NULL, log.p=TRUE, BPPARAM=SerialParam())
# This fits a linear model to each gene and performs a t-test for
# differential expression between groups.
{
    full.design <- model.matrix(~0 + groups)
    clust.vals <- levels(groups)
    colnames(full.design) <- clust.vals

    if (nrow(design)!=ncol(x)) {
        stop("'nrow(design)' is not equal to 'ncol(x)'")
    }
    if (any(is.intercept <- colnames(design) == "(Intercept)")) {
        design <- design[,!is.intercept,drop=FALSE]
        warning("automatically removed intercept column")
    }

    gene.names <- .setup_gene_names(gene.names, x, subset.row)

    # Pruning out NA 'groups'. Note that 'model.matrix'
    # automatically drops NA entries for 'full.design'.
    restricted <- !is.na(groups)
    if (!all(restricted)) {
        x <- x[,restricted,drop=FALSE]
        design <- design[restricted,,drop=FALSE]
    }

    # Linear dependencies will trigger errors in .ranksafe_QR.
    full.design <- cbind(full.design, design) 
    stats <- fitLinearModel(x, full.design, subset.row=subset.row, BPPARAM=BPPARAM)

    resid.df <- stats$residual.df
    if (resid.df <= 0L) {
        stop("no residual d.f. in design matrix for variance estimation")
    }

    coefficients <- stats$coefficients
    sigma2 <- stats$variance
    sigma2 <- pmax(sigma2, 1e-8) # avoid unlikely but possible problems with discreteness.

    collected.stats <- collected.pairs <- list()
    counter <- 1L

    # Running through every pair of groups.
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
        ref.coef <- coefficients[,h]
        for (tdex in seq_len(h-1L)) {
            target <- clust.vals[tdex]
            cur.lfc <- ref.coef - coefficients[,tdex]

            test.out <- .run_t_test(cur.lfc, lfit2$stdev.unscaled[tdex]^2*sigma2, resid.df, thresh.lfc=lfc)
            hvt.p <- .choose_leftright_pvalues(test.out$left, test.out$right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(test.out$right, test.out$left, direction=direction)

            if (std.lfc) {
                # Computing Cohen's D.
                is.zero <- cur.lfc == 0
                cur.lfc <- cur.lfc / sqrt(sigma2)
                cur.lfc[is.zero] <- 0
            }

            collected.stats[[counter]] <- list(
                .create_full_stats(cur.lfc, p=hvt.p, gene.names=gene.names, log.p=log.p),
                .create_full_stats(-cur.lfc, p=tvh.p, gene.names=gene.names, log.p=log.p)
            )
            collected.pairs[[counter]] <- DataFrame(first=c(host, target), second=c(target, host))
            counter <- counter + 1L
        }
    }

    output <- list(
        statistics=unlist(collected.stats, recursive=FALSE), 
        pairs=do.call(rbind, collected.pairs)
    )
    .reorder_pairwise_output(output, clust.vals)
}

###########################################################
# Internal functions (t-test calculations)
###########################################################

#' @importFrom stats pt
.run_t_test <- function(cur.lfc, cur.err, cur.df, thresh.lfc=0) 
# This runs the t-test given the relevant statistics, regardless of how
# they were computed (i.e., within blocks, or in a linear model).
{
    thresh.lfc <- abs(thresh.lfc)
    if (thresh.lfc==0) {
        cur.t <- cur.lfc/sqrt(cur.err)
        left <- pt(cur.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
        right <- pt(cur.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
    } else {
        # For one-sided tests, testing against the lfc threshold in the specified direction.
        # Note that if direction='up', only 'right' is used; nonetheless, 'left' is still calculated
        # using 'lower.t' so as to allow quick calculation of the p-value for the reversed contrast,
        # by simply swapping 'left' and 'right'. The same applies when direction='down'.
        lower.t <- (cur.lfc + thresh.lfc)/sqrt(cur.err)
        left <- pt(lower.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)

        upper.t <- (cur.lfc - thresh.lfc)/sqrt(cur.err)
        right <- pt(upper.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
    }

    list(left=left, right=right)
}

.add_log_values <- function(x, y)
# This performs log(exp(x) + exp(y)) in a numerically
# stable fashion, i.e., without actually running the 'exp's.
{
    pmax(x, y) + log1p(exp(-abs(x-y)))
}
