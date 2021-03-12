#' Perform pairwise Wilcoxon rank sum tests
#' 
#' Perform pairwise Wilcoxon rank sum tests between groups of cells, possibly after blocking on uninteresting factors of variation.
#' 
#' @param x A numeric matrix-like object of normalized (and possibly log-transformed) expression values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param direction A string specifying the direction of differences to be considered in the alternative hypothesis.
#' @param lfc Numeric scalar specifying the minimum log-fold change for one observation to be considered to be \dQuote{greater} than another.
#' @inheritParams pairwiseTTests
#' 
#' @details
#' This function performs Wilcoxon rank sum tests to identify differentially expressed genes (DEGs) between pairs of groups of cells.
#' A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of groups.
#' This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.
#' 
#' The effect size for each gene in each comparison is reported as the area under the curve (AUC).
#' Consider the distribution of expression values for gene X within each of two groups A and B.
#' The AUC is the probability that a randomly selected cell in A has a greater expression of X than a randomly selected cell in B. 
#' (Ties are assumed to be randomly broken.)
#' Concordance probabilities near 0 indicate that most observations in A are lower than most observations in B;
#' conversely, probabilities near 1 indicate that most observations in A are higher than those in B.
#' The Wilcoxon rank sum test effectively tests for significant deviations from an AUC of 0.5.
#' 
#' Wilcoxon rank sum tests are more robust to outliers and insensitive to non-normality, in contrast to t-tests in \code{\link{pairwiseTTests}}.
#' However, they take longer to run, the effect sizes are less interpretable, and there are more subtle violations of its assumptions in real data.
#' For example, the i.i.d. assumptions are unlikely to hold after scaling normalization due to differences in variance.
#' Also note that we approximate the distribution of the Wilcoxon rank sum statistic to deal with large numbers of cells and ties.
#' 
#' If \code{restrict} is specified, comparisons are only performed between pairs of groups in \code{restrict}.
#' This can be used to focus on DEGs distinguishing between a subset of the groups (e.g., closely related cell subtypes).
#'
#' If \code{exclude} is specified, comparisons are not performed between groups in \code{exclude}.
#' Similarly, if any entries of \code{groups} are \code{NA}, the corresponding cells are are ignored.
#' 
#' @section Direction and magnitude of the effect:
#' If \code{direction="any"}, two-sided Wilcoxon rank sum tests will be performed for each pairwise comparisons between groups of cells.
#' Otherwise, one-sided tests in the specified direction will be used instead.
#' This can be used to focus on genes that are upregulated in each group of interest, which is often easier to interpret.
#' 
#' To interpret the setting of \code{direction}, consider the DataFrame for group X, in which we are comparing to another group Y.
#' If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in group X compared to Y.
#' If \code{direction="down"}, genes will only be significant if they are downregulated in group X compared to Y.
#' See \code{?\link{wilcox.test}} for more details on the interpretation of one-sided Wilcoxon rank sum tests.
#'
#' Users can also specify a log-fold change threshold in \code{lfc} to focus on genes that exhibit large shifts in location.
#' This is equivalent to specifying the \code{mu} parameter in \code{\link{wilcox.test}} with some additional subtleties depending on \code{direction}:
#' \itemize{
#' \item If \code{direction="any"}, the null hypothesis is that the true shift lies in \code{[-lfc, lfc]}.
#' Two one-sided p-values are computed against the boundaries of this interval by shifting X's expression values to either side by \code{lfc},
#' and these are combined to obtain a (conservative) two-sided p-value.
#' \item If \code{direction="up"}, the null hypothesis is that the true shift is less than or equal to \code{lfc}.
#' A one-sided p-value is computed against the boundary of this interval.
#' \item If \code{direction="down"}, the null hypothesis is that the true shift is greater than or equal to \code{-lfc}.
#' A one-sided p-value is computed against the boundary of this interval.
#' }
#'
#' The AUC is conveniently derived from the U-statistic used in the test, which ensures that it is always consistent with the reported p-value.
#' An interesting side-effect is that the reported AUC is dependent on both the specified \code{lfc} and \code{direction}.
#' \itemize{
#' \item If \code{direction="up"}, the AUC is computed after shifting X's expression values to the left by the specified \code{lfc}.
#' An AUC above 0.5 means that X's values are \dQuote{greater} than Y's, even after shifting down the former by \code{lfc}.
#' This is helpful as a large AUC tells us that X and Y are well-separated by at least \code{lfc}.
#' However, an AUC below 0.5 cannot be interpreted as \dQuote{X is lower than Y}, only that \dQuote{X - lfc is lower than Y}.
#' \item If \code{direction="down"}, the AUC is computed after shifting X's expression values to the right by the specified \code{lfc}.
#' An AUC below 0.5 means that X's values are \dQuote{lower} than Y's, even after shifting up the former by \code{lfc}.
#' This is helpful as a low AUC tells us that X and Y are well-separated by at least \code{lfc}.
#' However, an AUC above 0.5 cannot be interpreted as \dQuote{X is greater than Y}, only that \dQuote{X + lfc is greater than Y}.
#' \item If \code{direction="any"}, the AUC is computed by averaging the AUCs obtained in each of the two one-sided tests, i.e., after shifting in each direction.
#' This considers an observation of Y to be tied with an observation of X if their absolute difference is less than \code{lfc}.
#' (Technically, the test procedure should also consider these to be ties to be fully consistent, but we have not done so for simplicity.)
#' The AUC can be interpreted as it would be for \code{lfc=0}, i.e., above 0.5 means that X is greater than Y and below 0.5 means that X is less than Y. 
#' }
#'
#' @section Blocking on uninteresting factors:
#' If \code{block} is specified, Wilcoxon tests are performed between groups of cells within each level of \code{block}.
#' For each pair of groups, the p-values for each gene across all levels of \code{block} are combined using Stouffer's Z-score method.
#' The reported AUC is also a weighted average of the AUCs across all levels.
#' 
#' The weight for a particular level of \code{block} is defined as \eqn{N_xN_y},
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
#' including the AUCs, p-values and false discovery rates.
#' 
#' The \code{pairs} element is a DataFrame with one row corresponding to each entry of \code{statistics}.
#' This contains the fields \code{first} and \code{second}, 
#' specifying the two groups under comparison in the corresponding DataFrame in \code{statistics}.
#' 
#' In each DataFrame in \code{statistics}, the AUC represents the probability of sampling a value in the \code{first} group greater than a random value from the \code{second} group.
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
#' @seealso
#' \code{\link{wilcox.test}}, on which this function is based.
#' 
#' \code{\link{combineMarkers}}, to combine pairwise comparisons into a single DataFrame per group.
#'
#' \code{\link{getTopMarkers}}, to obtain the top markers from each pairwise comparison.
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(sce)), centers=2) 
#' 
#' # Vanilla application:
#' out <- pairwiseWilcox(logcounts(sce), groups=kout$cluster)
#' out
#' 
#' # Directional and with a minimum log-fold change:
#' out <- pairwiseWilcox(logcounts(sce), groups=kout$cluster, 
#'     direction="up", lfc=0.2)
#' out
#' 
#' @export
#' @name pairwiseWilcox
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle .subset2index
.pairwiseWilcox <- function(x, groups, block=NULL, restrict=NULL, exclude=NULL, direction=c("any", "up", "down"),
    lfc=0, log.p=FALSE, gene.names=NULL, subset.row=NULL, BPPARAM=SerialParam())
{
    groups <- .setup_groups(groups, x, restrict=restrict, exclude=exclude) 
    direction <- match.arg(direction)

    # Actual calculations occur inside another function, for symmetry with pairwiseTTests.
    .blocked_wilcox(x, subset.row, groups, block=block, direction=direction, 
        lfc=lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
}

#' @export
#' @rdname pairwiseWilcox
setGeneric("pairwiseWilcox", function(x, ...) standardGeneric("pairwiseWilcox"))

#' @export
#' @rdname pairwiseWilcox
setMethod("pairwiseWilcox", "ANY", .pairwiseWilcox)

#' @export
#' @rdname pairwiseWilcox
#' @importFrom SummarizedExperiment assay
setMethod("pairwiseWilcox", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .pairwiseWilcox(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname pairwiseWilcox
#' @importFrom SingleCellExperiment colLabels
setMethod("pairwiseWilcox", "SingleCellExperiment", function(x, groups=colLabels(x, onAbsence="error"), ...) {
    callNextMethod(x=x, groups=groups, ...)
})

###########################################################
# Internal functions (blocking)
###########################################################

#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom stats pnorm 
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom beachmat rowBlockApply
.blocked_wilcox <- function(x, subset.row, groups, block=NULL, direction="any", gene.names=NULL, 
    lfc=0, log.p=TRUE, BPPARAM=SerialParam())
{
    if (is.null(block)) {
        block <- list(`1`=seq_len(ncol(x)))
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- split(seq_along(block), block)
    }

    gene.names <- .setup_gene_names(gene.names, x, subset.row)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    # Setting up the parallelization strategy.
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Computing across blocks.
    group.vals <- levels(groups)
    nblocks <- length(block)
    all.stats <- all.ties <- all.n <- vector("list", nblocks)

    for (b in seq_along(block)) {
        chosen <- block[[b]]
        cur.groups <- groups[chosen]
        all.n[[b]] <- as.vector(table(cur.groups))
        names(all.n[[b]]) <- group.vals
        
        by.group <- split(chosen - 1L, cur.groups)
        bpl.out <- rowBlockApply(x, FUN=overlap_exprs, groups=by.group, lfc=lfc, BPPARAM=BPPARAM)
        raw.stats <- lapply(bpl.out, "[[", i=1)
        raw.ties <- lapply(bpl.out, "[[", i=2)

        cons.stats <- cons.ties <- vector("list", length(group.vals))
        names(cons.stats) <- names(cons.ties) <- group.vals
        for (i in seq_along(cons.stats)) {
            cons.stats[[i]] <- do.call(rbind, lapply(raw.stats, "[[", i=i))
            cons.ties[[i]] <- do.call(rbind, lapply(raw.ties, "[[", i=i))
            ncols <- ncol(cons.stats[[i]])
            colnames(cons.stats[[i]]) <- colnames(cons.ties[[i]]) <- group.vals[seq_len(ncols)]
        }

        all.stats[[b]] <- cons.stats
        all.ties[[b]] <- cons.ties
    }

    if (lfc==0) {
        STATFUN <- .generate_nolfc_wilcox(all.n, all.stats, all.ties, direction)
    } else {
        STATFUN <- .generate_lfc_wilcox(all.n, all.stats, all.ties, direction)
    }

    # This looks at every level of the blocking factor and performs
    # Wilcoxon tests between pairs of groups within each blocking level.
    .pairwise_blocked_template(group.vals, nblocks=length(block), direction=direction, 
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, effect.name="AUC",
        BPPARAM=BPPARAM)
}

##########################
### Internal functions ###
##########################

#' @importFrom stats pnorm
.generate_nolfc_wilcox <- function(all.n, all.stats, all.ties, direction) {
    force(all.n)
    force(all.stats)
    force(all.ties)
    force(direction)

    function(b, host, target) {
        host.n <- as.double(all.n[[b]][[host]]) # numeric conversion to avoid overflow.
        target.n <- as.double(all.n[[b]][[target]])
        cur.prod <- host.n * target.n

        effect <- all.stats[[b]][[host]][,target]
        auc <- effect/cur.prod
        output <- list(forward=auc, reverse=1 - auc, weight=cur.prod)

        # 'cur.prod' is still nominally integer; we use 1.5 to avoid
        # numerical imprecision upon an exact comparison.
        output$valid <- cur.prod > 1.5

        # Approximate Wilcoxon with continuity correction: ripped straight from wilcox.test() in stats.
        z <- effect - cur.prod/2
        SIGMA <- .get_sigma(host.n, target.n, all.ties[[b]][[host]][,target])

        # Always dealing with one-sided tests, so we fix the correction at 0.5
        # rather than allowing it to be zero for direction='any' (as in wilcox.test()).
        CORRECTION <- 0.5
        output$left <- pnorm((z + CORRECTION)/SIGMA, log.p=TRUE)
        output$right <- pnorm((z - CORRECTION)/SIGMA, log.p=TRUE, lower.tail=FALSE)

        output
    }
}

#' @importFrom stats pnorm
.generate_lfc_wilcox <- function(all.n, all.stats, all.ties, direction) {
    force(all.n)
    force(all.stats)
    force(all.ties)
    force(direction)

    function(b, host, target) {
        host.n <- as.double(all.n[[b]][[host]]) # numeric conversion to avoid overflow.
        target.n <- as.double(all.n[[b]][[target]])
        cur.prod <- host.n * target.n

        # Minus, in that host's values have 'lfc' subtracted from them (i.e., null is +lfc).
        # Added, in that host's values have 'lfc' added to them (i.e., null is -lfc).
        minus.effect <- all.stats[[b]][[host]][,target]
        added.effect <- all.stats[[b]][[target]][,host]

        if (direction=="any") {
            # Taking the average to ensure that the AUC is interpretable around 0.5, see Details.
            effect <- minus.effect/2 + added.effect/2
            auc <- effect/cur.prod
            output <- list(forward=auc, reverse=1-auc)
        } else if (direction=="up") {
            output <- list(forward=minus.effect/cur.prod, reverse=1-added.effect/cur.prod)
        } else {
            output <- list(forward=added.effect/cur.prod, reverse=1-minus.effect/cur.prod)
        }
        output$weight <- cur.prod

        # 'cur.prod' is still nominally integer; we use 1.5 to avoid
        # numerical imprecision upon an exact comparison.
        output$valid <- cur.prod > 1.5

        minus.z <- minus.effect - cur.prod/2
        minus.SIGMA <- .get_sigma(host.n, target.n, all.ties[[b]][[host]][,target])
        added.z <- added.effect - cur.prod/2
        added.SIGMA <- .get_sigma(host.n, target.n, all.ties[[b]][[target]][,host])

        CORRECTION <- 0.5

        # For one-sided tests, testing against the lfc threshold in the specified direction.
        # Note that if direction='up', only 'right' is used; nonetheless, 'left' is still calculated
        # so as to allow quick calculation of the p-value for the reversed contrast,
        # by simply swapping 'left' and 'right'. The same applies when direction='down'.
        output$left <- pnorm((added.z + CORRECTION)/added.SIGMA, log.p=TRUE)
        output$right <- pnorm((minus.z - CORRECTION)/minus.SIGMA, log.p=TRUE, lower.tail=FALSE)

        output
    }
}

.get_sigma <- function(host.n, target.n, cur.ties) {
    s2 <- (host.n * target.n/12) * ((host.n + target.n + 1) - cur.ties/((host.n + target.n) * (host.n + target.n - 1)))
    pmax(sqrt(s2), 1e-8)
}
