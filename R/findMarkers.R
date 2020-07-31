#' Find marker genes
#' 
#' Find candidate marker genes for groups of cells (e.g., clusters) 
#' by testing for differential expression between pairs of groups.
#' 
#' @param x A numeric matrix-like object of expression values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#' This is expected to be normalized log-expression values for most tests - see Details.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param groups A vector of length equal to \code{ncol(x)}, 
#' specifying the group to which each cell is assigned.
#' If \code{x} is a \linkS4class{SingleCellExperiment}, this defaults to \code{\link{colLabels}(x)} if available.
#' @param test.type String specifying the type of pairwise test to perform -
#' a t-test with \code{"t"}, a Wilcoxon rank sum test with \code{"wilcox"}, 
#' or a binomial test with \code{"binom"}.
#' @inheritParams combineMarkers
#' @param log.p A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param row.data A \linkS4class{DataFrame} containing additional row metadata for each gene in \code{x},
#' to be included in each of the output DataFrames.
#' This should generally have row names identical to those of \code{x}.
#' 
#' Alternatively, a list containing one such DataFrame per level of \code{groups}, 
#' where each DataFrame contains group-specific metadata for each gene to be included in the appropriate output DataFrame.
#' @param add.summary Logical scalar indicating whether statistics from \code{\link{summaryMarkerStats}} should be added.
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
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param assay.type A string specifying which assay values to use, usually \code{"logcounts"}.
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
#' See \code{?\link{combineMarkers}} for more details on the output format.
#'
#' If \code{row.data} is provided, the additional fields are added to the front of the DataFrame for each cluster.
#' If \code{add.summary=TRUE}, extra statistics for each cluster are also computed and added.
#' 
#' Any log-fold changes are reported as differences in average \code{x} between groups
#' (usually in base 2, depending on the transformation applied to \code{x}).
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
#' \code{\link{summaryMarkerStats}}, to incorporate additional summary statistics per cluster.
#' 
#' \code{\link{getMarkerEffects}}, to easily extract a matrix of effect sizes from each DataFrame.
#'
#' @examples
#' library(scuttle)
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

#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom BiocGenerics cbind
#' @importFrom scuttle .bpNotSharedOrUp
.findMarkers <- function(x, groups, test.type=c("t", "wilcox", "binom"), ..., 
    pval.type=c("any", "some", "all"), min.prop=NULL, log.p=FALSE, full.stats=FALSE, 
    sorted=TRUE, row.data=NULL, add.summary=FALSE, BPPARAM=SerialParam())
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

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    fit <- FUN(x, groups, ..., log.p=TRUE, BPPARAM=BPPARAM)

    output <- combineMarkers(fit$statistics, fit$pairs, pval.type=pval.type, min.prop=min.prop, 
        log.p.in=TRUE, log.p.out=log.p, full.stats=full.stats, pval.field="log.p.value", 
        effect.field=effect.field, sorted=sorted, BPPARAM=BPPARAM)

    if (add.summary) {
        row.data <- summaryMarkerStats(x, groups, row.data=row.data, BPPARAM=BPPARAM)
    }

    .add_row_data(output, row.data, match.names=sorted)
}

#' @importClassesFrom S4Vectors DataFrame
.add_row_data <- function(output, row.data, match.names) {
    if (is.null(row.data)) {
        return(output)
    }

    for (i in names(output)) {
        current <- output[[i]]

        if (is.data.frame(row.data) || is(row.data, "DataFrame")) {
            rd <- row.data
        } else {
            if (!i %in% names(row.data)) {
                stop("list-like 'row.data' should be named with the levels of 'groups'")
            }
            rd <- row.data[[i]]
        }

        rn <- rownames(current)
        if (match.names) {
            if (is.null(rn) || !identical(sort(rn), sort(rownames(rd)))) {
                stop("inconsistent or NULL row names for 'row.data' and result tables")
            }
            rd <- rd[rn,,drop=FALSE]
        } else if (!identical(rn, rownames(rd))) {
            stop("inconsistent row names for 'row.data' and result tables")
        }

        output[[i]] <- cbind(rd, current)
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
setMethod("findMarkers", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .findMarkers(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname findMarkers
#' @importFrom SingleCellExperiment colLabels
setMethod("findMarkers", "SingleCellExperiment", function(x, groups=colLabels(x, onAbsence="error"), ...) { 
    callNextMethod(x=x, groups=groups, ...)
})
