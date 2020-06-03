#' Get top markers
#'
#' Obtain the top markers for each pairwise comparison between clusters, or for each cluster.
#'
#' @inheritParams combineMarkers
#' @param pval.field String specifying the column of each DataFrame in \code{de.lists} to use to identify top markers.
#' Smaller values are assigned higher rank.
#' @param n Integer scalar specifying the number of markers to obtain from each pairwise comparison, if \code{pairwise=FALSE}.
#'
#' Otherwise, the number of top genes to take from each cluster's combined marker set, see Details.
#' @param pairwise Logical scalar indicating whether top markers should be returned for every pairwise comparison.
#' If \code{FALSE}, one marker set is returned for every cluster.
#' @param pval.type String specifying how markers from pairwise comparisons are to be combined if \code{pairwise=FALSE}.
#' This has the same effect as \code{pval.type} in \code{\link{combineMarkers}}.
#' @param fdr.field String specifying the column containing the adjusted p-values.
#' @param fdr.threshold Numeric scalar specifying the FDR threshold for filtering.
#' If \code{NULL}, no filtering is performed on the FDR.
#' @param ... Further arguments to pass to \code{\link{combineMarkers}} if \code{pairwise=FALSE}.
#'
#' @return If \code{pairwise=TRUE}, a \linkS4class{List} of Lists of character vectors is returned.
#' Each element of the outer list corresponds to cluster X, each element of the inner list corresponds to another cluster Y,
#' and each character vector specifies the marker genes that distinguish X from Y.
#'
#' If \code{pairwise=FALSE}, a List of character vectors is returned.
#' Each character vector contains the marker genes that distinguish X from any, some or all other clusters,
#' depending on \code{combine.type}.
#'
#' @author Aaron Lun
#'
#' @details
#' This is a convenience utility that converts the results of pairwise comparisons into a marker list
#' that can be used in downstream functions, e.g., as the marker sets in \pkg{SingleR}.
#' By default, it returns a list of lists containing the top genes for every pairwise comparison,
#' which is useful for feature selection to select genes distinguishing between closely related clusters.
#' The top \code{n} genes are chosen with adjusted p-values below \code{fdr.threshold}.
#'
#' If \code{pairwise=FALSE}, \code{\link{combineMarkers}} is called on \code{de.lists} and \code{pairs}
#' to obtain a per-cluster ranking of genes from all pairwise comparisons involving that cluster.
#' If \code{pval.type="any"}, the top genes with \code{Top} values no greater than \code{n} are retained; 
#' this is equivalent to taking the union of the top \code{n} genes from each pairwise comparison for each cluster.
#' Otherwise, the top \code{n} genes with the smallest p-values are retained.
#' In both cases, genes are further filtered by \code{fdr.threshold}.
#' 
#' @seealso
#' \code{\link{pairwiseTTests}} and friends, to obtain \code{de.lists} and \code{pairs}.
#'
#' \code{\link{combineMarkers}}, for another function that consolidates pairwise DE comparisons.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(sce)), centers=3) 
#' 
#' out <- pairwiseTTests(logcounts(sce), 
#'      groups=paste0("Cluster", kout$cluster))
#'
#' # Getting top pairwise markers:
#' top <- getTopMarkers(out$statistics, out$pairs)
#' top[[1]]
#' top[[1]][[2]]
#' 
#' # Getting top per-cluster markers:
#' top <- getTopMarkers(out$statistics, out$pairs, pairwise=FALSE)
#' top[[1]]
#' @export
#' @importFrom S4Vectors List
#' @importFrom utils head
getTopMarkers <- function(de.lists, pairs, n=10, pval.field="p.value", fdr.field="FDR", 
    pairwise=TRUE, pval.type=c("any", "some", "all"), fdr.threshold=0.05, ...)
{
    markers <- List()
    all.labels <- sort(unique(c(pairs$first, pairs$second)))

    .sigrows <- function(stats) {
        if (!is.null(fdr.threshold)) {
            cur.fdr <- stats[,fdr.field]
            stats <- stats[cur.fdr <= fdr.threshold & !is.na(cur.fdr),,drop=FALSE]
        }
        stats
    }

    if (pairwise) {
        for (first in all.labels) {
            cur.markers <- List()
            for (second in all.labels) {
                chosen <- which(pairs$first==first & pairs$second==second)

                if (!length(chosen)) {
                    cur.markers[[second]] <- character(0)
                } else if (length(chosen)!=1L){ 
                    stop(sprintf("multiple entries in 'pairs' for '%s' vs '%s'", first, second))
                } else {
                    cur.stats <- .sigrows(de.lists[[chosen]])
                    o <- order(cur.stats[[pval.field]])
                    cur.markers[[second]] <- rownames(cur.stats)[head(o, n)]
                }
            }

            markers[[first]] <- cur.markers
        }
    } else {
        pval.type <- match.arg(pval.type)
        combined <- combineMarkers(de.lists, pairs, ..., pval.field=pval.field,
            effect.field=NULL, pval.type=pval.type, sorted=TRUE)

        for (i in names(combined)) {
            cur.stats <- .sigrows(combined[[i]])
            if (pval.type=="any") {
                markers[[i]] <- rownames(cur.stats)[cur.stats$Top <= n]
            } else {
                markers[[i]] <- head(rownames(cur.stats), n)
            }
        }
    }

    markers
}
