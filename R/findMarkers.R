#' Find marker genes
#' 
#' Find candidate marker genes for groups of cells (e.g., clusters) by testing for differential expression between pairs of groups.
#' 
#' @param x A numeric matrix-like object of expression values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#' This is expected to be normalized log-expression values for most tests - see Details.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param groups A vector of group assignments for all cells.
#' @param test.type String specifying the type of pairwise test to perform -
#' a t-test with \code{"t"}, a Wilcoxon rank sum test with \code{"wilcox"}, 
#' or a binomial test with \code{"binom"}.
#' @inheritParams combineMarkers
#' @param log.p A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param row.data A \linkS4class{DataFrame} containing additional row metadata for each gene in \code{x},
#' to be included in each of the output DataFrames.
#' If \code{sorted=TRUE}, this should have the same row names as the output of \code{\link{combineMarkers}}
#' (usually \code{rownames(x)}, see the \code{gene.names} argument in functions like \code{\link{pairwiseTTests}}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY method:
#' \itemize{
#' \item For \code{test.type="t"}, further arguments to pass to \code{\link{pairwiseTTests}}.
#' \item For \code{test.type="wilcox"}, further arguments to pass to \code{\link{pairwiseWilcox}}.
#' \item For \code{test.type="binom"}, further arguments to pass to \code{\link{pairwiseBinom}}.
#' }
#' Common arguments for all testing functions include \code{gene.names}, \code{direction}, 
#' \code{block} and \code{BPPARAM}.
#' Test-specific arguments are also supported for the appropriate \code{test.type}.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use, usually \code{"logcounts"}.
#' @param get.spikes See \code{?"\link{scran-gene-selection}"}.
#' 
#' @details
#' This function provides a convenience wrapper for marker gene identification between groups of cells,
#' based on running \code{\link{pairwiseTTests}} or related functions and passing the result to \code{\link{combineMarkers}}.
#' All of the arguments above are supplied directly to one of these two functions -
#' refer to the relevant function's documentation for more details.
#' 
#' If \code{x} contains log-normalized expression values generated with a pseudo-count of 1,
#' it can be used in any of the pairwise testing procedures.
#' If \code{x} is scale-normalized but not log-transformed, it can be used with \code{test.type="wilcox"} and \code{test.type="binom"}.
#' If \code{x} contains raw counts, it can only be used with \code{test.type="binom"}.
#' 
#' Note that \code{log.p} only affects the combined p-values and FDRs.
#' If \code{full.stats=TRUE}, the p-values for each individual pairwise comparison will always be log-transformed,
#' regardless of the value of \code{log.p}.
#' Log-transformed p-values and FDRs are reported using the natural base.
#'
#' The choice of \code{pval.type} determines whether the highly ranked genes are those that are DE between the current group and:
#' \itemize{
#' \item any other group (\code{"any"})
#' \item all other groups (\code{"all"})
#' \item some other groups (\code{"some"})
#' }
#' See \code{?\link{combineMarkers}} for more details.
#'
#' @return 
#' A named list of \linkS4class{DataFrame}s, each of which contains a sorted marker gene list for the corresponding group.
#' In each DataFrame, the top genes are chosen to enable separation of that group from all other groups.
#' Log-fold changes are reported as differences in average \code{x} between groups
#' (usually in base 2, depending on the transformation applied to \code{x}).
#' 
#' See \code{?\link{combineMarkers}} for more details on the output format.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{pairwiseTTests}},
#' \code{\link{pairwiseWilcox}},
#' \code{\link{pairwiseBinom}},
#' for the underlying functions that compute the pairwise DE statistics.
#'
#' \code{\link{combineMarkers}}, to combine pairwise statistics into a single marker list per cluster.
#' 
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay, only using k-means for convenience.
#' kout <- kmeans(t(logcounts(sce)), centers=4) 
#'
#' out <- findMarkers(sce, groups=kout$cluster)
#' names(out)
#' out[[1]]
#'
#' # More customization of the tests:
#' out <- findMarkers(sce, groups=kout$cluster, test.type="wilcox")
#' out[[1]]
#'
#' out <- findMarkers(sce, groups=kout$cluster, lfc=1, direction="up")
#' out[[1]]
#'
#' out <- findMarkers(sce, groups=kout$cluster, pval.type="all")
#' out[[1]]
#'
#' @name findMarkers
NULL

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocGenerics cbind
.findMarkers <- function(x, groups, test.type=c("t", "wilcox", "binom"), ..., 
    pval.type=c("any", "some", "all"), min.prop=NULL, log.p=FALSE, 
    full.stats=FALSE, sorted=TRUE, row.data=NULL) 
{
    test.type <- match.arg(test.type)
    if (test.type=="t") {
        FUN <- pairwiseTTests
        effect.field <- "logFC"
    } else if (test.type=="wilcox") {
        FUN <- pairwiseWilcox
        effect.field <- "AUC"
    } else {
        FUN <- pairwiseBinom
        effect.field <- "logFC"
    }

    fit <- FUN(x, groups, ..., log.p=TRUE)
    output <- combineMarkers(fit$statistics, fit$pairs, pval.type=pval.type, min.prop=min.prop, 
        log.p.in=TRUE, log.p.out=log.p, full.stats=full.stats, pval.field="log.p.value", 
        effect.field=effect.field, sorted=sorted)

    if (!is.null(row.data)) {
        for (i in seq_along(output)) {
            if (sorted) {
                if (!identical(sort(rownames(output[[i]])), sort(rownames(row.data)))) {
                    stop("non-identical row names for 'row.data' and result tables")
                }
                cur.row.data <- row.data[rownames(output[[i]]),,drop=FALSE]
            } else {
                cur.row.data <- row.data
                rownames(cur.row.data) <- rownames(output[[i]])
            }
            output[[i]] <- cbind(cur.row.data, output[[i]])
        }
    }
    output
}

#' @export
#' @rdname findMarkers
setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

#' @export
#' @rdname findMarkers
setMethod("findMarkers", "ANY", .findMarkers)

#' @export
#' @rdname findMarkers
#' @importFrom SummarizedExperiment assay
setMethod("findMarkers", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .findMarkers(assay(x, i=assay.type), ..., subset.row=subset.row)
})
