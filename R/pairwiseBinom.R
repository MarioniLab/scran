#' Perform pairwise binomial tests
#'
#' Perform pairwise binomial tests between groups of cells, 
#' possibly after blocking on uninteresting factors of variation.
#' 
#' @param x A numeric matrix-like object of counts,
#' where each column corresponds to a cell and each row corresponds to a gene.
#' @param direction A string specifying the direction of effects to be considered for the alternative hypothesis.
#' @param lfc Numeric scalar specifying the minimum absolute log-ratio in the proportion of expressing genes between groups.
#' @param threshold Numeric scalar specifying the value below which a gene is presumed to be not expressed.
#' @inheritParams pairwiseTTests
#' 
#' @details
#' This function performs exact binomial tests to identify marker genes between pairs of groups of cells.
#' Here, the null hypothesis is that the proportion of cells expressing a gene is the same between groups.
#' A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of groups.
#' This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.
#' 
#' Effect sizes for each comparison are reported as log2-fold changes in the proportion of expressing cells in one group over the proportion in another group.
#' We add a pseudo-count that squeezes the log-FCs towards zero to avoid undefined values when one proportion is zero.
#' This is closely related to but somewhat more interpretable than the log-odds ratio,
#' which would otherwise be the more natural statistic for a proportion-based test.
#'
#' If \code{restrict} is specified, comparisons are only performed between pairs of groups in \code{restrict}.
#' This can be used to focus on DEGs distinguishing between a subset of the groups (e.g., closely related cell subtypes).
#' Similarly, if any entries of \code{groups} are \code{NA}, the corresponding cells are are ignored.
#' 
#' \code{x} can be a count matrix or any transformed counterpart where zeroes remain zero and non-zeroes remain non-zero.
#' This is true of any scaling normalization and monotonic transformation like the log-transform.
#' If the transformation breaks this rule, some adjustment of \code{threshold} is necessary.
#' 
#' A consequence of the transformation-agnostic behaviour of this function is that it will not respond to normalization.
#' Differences in library size will not be considered by this function.
#' However, this is not necessarily problematic for marker gene detection -
#' users can treat this as \emph{retaining} information about the total RNA content, analogous to spike-in normalization.
#' 
#' @section Direction and magnitude of the effect:
#' If \code{direction="any"}, two-sided binomial tests will be performed for each pairwise comparisons between groups of cells.
#' For other \code{direction}, one-sided tests in the specified direction will be used instead. 
#' This can be used to focus on genes that are upregulated in each group of interest, which is often easier to interpret.
#' 
#' In practice, the two-sided test is approximated by combining two one-sided tests using a Bonferroni correction.
#' This is done for various logistical purposes;
#' it is also the only way to combine p-values across blocks in a direction-aware manner.
#' As a result, the two-sided p-value reported here will not be the same as that from \code{\link{binom.test}}.
#' In practice, they are usually similar enough that this is not a major concern.
#' 
#' To interpret the setting of \code{direction}, consider the DataFrame for group X, in which we are comparing to another group Y.
#' If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in group X compared to Y.
#' If \code{direction="down"}, genes will only be significant if they are downregulated in group X compared to Y.
#' See \code{?\link{binom.test}} for more details on the interpretation of one-sided Wilcoxon rank sum tests.
#'
#' The magnitude of the log-fold change in the proportion of expressing cells can also be tested by setting \code{lfc}.
#' By default, \code{lfc=0} meaning that we will reject the null upon detecting any difference in proportions.
#' If this is set to some other positive value, the null hypothesis will change depending on \code{direction}:
#' \itemize{
#' \item If \code{direction="any"}, the null hypothesis is that the true log-fold change in proportions lies within \code{[-lfc, lfc]}.
#' To be conservative, we perform one-sided tests against the boundaries of this interval, and combine them to obtain a two-sided p-value.
#' \item If \code{direction="up"}, the null hypothesis is that the true log-fold change is less than \code{lfc}.
#' A one-sided p-value is computed against the boundary of this interval.
#' \item If \code{direction="down"}, the null hypothesis is that the true log-fold change is greater than \code{-lfc}.
#' A one-sided p-value is computed against the boundary of this interval.
#' }
#' 
#' @section Blocking on uninteresting factors:
#' If \code{block} is specified, binomial tests are performed between groups of cells within each level of \code{block}.
#' For each pair of groups, the p-values for each gene across 
#' all levels of \code{block} are combined using Stouffer's weighted Z-score method.
#' 
#' The weight for the p-value in a particular level of \code{block} is defined as \eqn{N_x + N_y},
#' where \eqn{N_x} and \eqn{N_y} are the number of cells in groups X and Y, respectively, for that level. 
#' This means that p-values from blocks with more cells will have a greater contribution to the combined p-value for each gene.
#' 
#' When combining across batches, one-sided p-values in the same direction are combined first.
#' Then, if \code{direction="any"}, the two combined p-values from both directions are combined.
#' This ensures that a gene only receives a low overall p-value if it changes in the same direction across batches.
#'
#' When comparing two groups, blocking levels are ignored if no p-value was reported, e.g., if there were insufficient cells for a group in a particular level. 
#' If all levels are ignored in this manner, the entire comparison will only contain \code{NA} p-values and a warning will be emitted.
#' 
#' @return
#' A list is returned containing \code{statistics} and \code{pairs}.
#' 
#' The \code{statistics} element is itself a list of \linkS4class{DataFrame}s.
#' Each DataFrame contains the statistics for a comparison between a pair of groups,
#' including the overlap proportions, p-values and false discovery rates.
#' 
#' The \code{pairs} element is a DataFrame with one row corresponding to each entry of \code{statistics}.
#' This contains the fields \code{first} and \code{second}, 
#' specifying the two groups under comparison in the corresponding DataFrame in \code{statistics}.
#' 
#' In each DataFrame in \code{statistics}, the log-fold change represents the log-ratio of the proportion of expressing cells in the \code{first} group compared to the expressing proportion in the \code{second} group.
#' 
#' @author
#' Aaron Lun
#' 
#' @references
#' Whitlock MC (2005). 
#' Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
#' \emph{J. Evol. Biol.} 18, 5:1368-73.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(sce)), centers=2) 
#' 
#' # Vanilla application:
#' out <- pairwiseBinom(logcounts(sce), groups=kout$cluster)
#' out
#' 
#' # Directional and with a minimum log-fold change:
#' out <- pairwiseBinom(logcounts(sce), groups=kout$cluster, 
#'     direction="up", lfc=1)
#' out
#'
#' @seealso
#' \code{\link{binom.test}} and \code{\link{binomTest}}, on which this function is based.
#'
#' \code{\link{combineMarkers}}, to combine pairwise comparisons into a single DataFrame per group.
#'
#' \code{\link{getTopMarkers}}, to obtain the top markers from each pairwise comparison.
#' @export
#' @name pairwiseBinom
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle .subset2index
.pairwiseBinom <- function(x, groups, block=NULL, restrict=NULL, exclude=NULL, direction=c("any", "up", "down"),
    threshold=1e-8, lfc=0, log.p=FALSE, gene.names=NULL, 
    subset.row=NULL, BPPARAM=SerialParam())
{
    groups <- .setup_groups(groups, x, restrict=restrict, exclude=exclude)
    direction <- match.arg(direction)

    # Actual calculations occur inside another function, for symmetry with pairwiseTTests.
    .blocked_binom(x, subset.row, groups, block=block, direction=direction, 
        gene.names=gene.names, log.p=log.p, threshold=threshold, lfc=lfc, BPPARAM=BPPARAM)
}

#' @export
#' @rdname pairwiseBinom
setGeneric("pairwiseBinom", function(x, ...) standardGeneric("pairwiseBinom"))

#' @export
#' @rdname pairwiseBinom
setMethod("pairwiseBinom", "ANY", .pairwiseBinom)

#' @export
#' @rdname pairwiseBinom
#' @importFrom SummarizedExperiment assay
setMethod("pairwiseBinom", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .pairwiseBinom(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname pairwiseBinom
#' @importFrom SingleCellExperiment colLabels
setMethod("pairwiseBinom", "SingleCellExperiment", function(x, groups=colLabels(x, onAbsence="error"), ...) {
    callNextMethod(x=x, groups=groups, ...)
})

###########################################################
# Internal functions (blocking)
###########################################################

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom stats pbinom
#' @importFrom scuttle .bpNotSharedOrUp numDetectedAcrossCells
#' @importFrom SummarizedExperiment assay
.blocked_binom <- function(x, subset.row, groups, block=NULL, direction="any", gene.names=NULL, log.p=TRUE, 
	threshold=1e-8, lfc=0, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# binomial tests between pairs of groups within each blocking level.
{
    if (is.null(block)) {
        block <- list(`1`=seq_len(ncol(x)))
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- split(seq_along(block), block)
    }

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Computing across blocks.
    group.vals <- levels(groups)
    nblocks <- length(block)
    all.nzero <- all.n <- vector("list", nblocks)

    for (b in seq_along(block)) {
        chosen <- block[[b]]
        cur.groups <- groups[chosen]
        all.n[[b]] <- as.vector(table(cur.groups))
        names(all.n[[b]]) <- group.vals

        raw.nzero <- numDetectedAcrossCells(x[,chosen,drop=FALSE], subset.row=subset.row, 
            ids=cur.groups, detection_limit=threshold, BPPARAM=BPPARAM)
        raw.nzero <- assay(raw.nzero)

        if (any(!group.vals %in% colnames(raw.nzero))) {
            # Handle missing levels gracefully.
            tmp <- matrix(0L, nrow=nrow(raw.nzero), ncol=length(group.vals),
                dimnames=list(rownames(raw.nzero), group.vals))
            tmp[,colnames(raw.nzero)] <- raw.nzero
            raw.nzero <- tmp
        }

        all.nzero[[b]] <- raw.nzero
    }

    if (lfc==0) {
        STATFUN <- .generate_nolfc_binom(all.n, all.nzero)
    } else {
        STATFUN <- .generate_lfc_binom(all.n, all.nzero, direction, lfc)
    }

    gene.names <- .setup_gene_names(gene.names, x, subset.row)

    .pairwise_blocked_template(group.vals, nblocks=length(block), direction=direction, 
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, effect.name="logFC",
        BPPARAM=BPPARAM)
}

##########################
### Internal functions ###
##########################

#' @importFrom stats pbinom
.generate_nolfc_binom <- function(all.n, all.nzero) {
    force(all.n)
    force(all.nzero)

    function(b, host, target) {    
        host.nzero <- all.nzero[[b]][,host]
        target.nzero <- all.nzero[[b]][,target]
        host.n <- all.n[[b]][[host]]
        target.n <- all.n[[b]][[target]]

        p <- host.n/(host.n + target.n)
        size <- host.nzero + target.nzero
        effect <- .compute_binom_effect(host.nzero, host.n, target.nzero, target.n)

        # Not exactly equal to binom.test(); the two-sided p-value from binom.test()
        # cannot be performed by any combination of the one-sided p-values. This 
        # makes it impossible to behave with directional Stouffer's method in 
        # .pairwise_blocking_template(), and just generally gums up the works.
        list(
            forward=effect, 
            reverse=-effect,
            weight=as.double(host.n) + as.double(target.n),
            valid=host.n > 0L && target.n > 0L,
            left=pbinom(host.nzero, size, p, log.p=TRUE),
            right=pbinom(host.nzero - 1, size, p, lower.tail=FALSE, log.p=TRUE)
        )
    }
}

# Log-fold change in proportions, mimic edgeR::cpm().
.compute_binom_effect <- function(host.nzero, host.n, target.nzero, target.n) {
    mean.lib <- mean(c(host.n, target.n))
    pseudo.host <- 1 * host.n/mean.lib
    pseudo.target <- 1 * target.n/mean.lib
    unname(log2((host.nzero + pseudo.host)/(host.n + 2 * pseudo.host))
        - log2((target.nzero + pseudo.target)/(target.n + 2 * pseudo.target)))
}

#' @importFrom stats pbinom
.generate_lfc_binom <- function(all.n, all.nzero, direction, lfc) {
    force(all.n)
    force(all.nzero)
    force(direction)
    fold <- 2^lfc

    function(b, host, target) {    
        host.nzero <- all.nzero[[b]][,host]
        target.nzero <- all.nzero[[b]][,target]
        host.n <- all.n[[b]][[host]]
        target.n <- all.n[[b]][[target]]

        # Converting log-fold change in the **proportions** into **probabilities**.
        # Let p_1 and p_2 be the probability of non-zero in group 1 and 2 respectively.
        # Let's say that p_1 = fold * p_2 for some fold > 1.
        # Given a non-zero value, the probability that it comes from group 1 is 
        # (p_1 * n_1) / (p_1 * n_1 + p_2 * n_2), which collapses to `p.right`
        # (i.e., probability of more non-zeros, hence the right side of the distribution).
        # Calculation of `p.left` follows the opposite premise that p_1 = p_2 / fold, 
        # i.e., the other side of the composite null hypothesis.
        p.left <- host.n/fold / (target.n + host.n/fold)
        p.right <- host.n*fold / (target.n + host.n*fold)
        size <- host.nzero + target.nzero
        effect <- .compute_binom_effect(host.nzero, host.n, target.nzero, target.n)

        output <- list(
            forward=effect, 
            reverse=-effect,
            weight=as.double(host.n) + as.double(target.n),
            valid=host.n > 0L && target.n > 0L
        )

        output$left <- pbinom(host.nzero, size, p.left, log.p=TRUE)
        output$right <- pbinom(host.nzero - 1, size, p.right, lower.tail=FALSE, log.p=TRUE)

        output
    }
}
