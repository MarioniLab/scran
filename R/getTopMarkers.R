#' Get top markers
#'
#' Obtain the top markers for each pairwise comparison between clusters, or for each cluster.
#'
#' @inheritParams combineMarkers
#' @param field String specifying the column of each DataFrame in \code{de.lists} to use to identify top markers.
#' Smaller values are assigned higher rank.
#' @param n Integer scalar specifying the number of markers to obtain from each pairwise comparison.
#' @param pairwise Logical scalar indicating whether top markers should be returned for every pairwise comparison.
#' If \code{FALSE}, one marker set is returned for every cluster.
#' @param combine.type String specifying how markers from pairwise comparisons are to be combined if \code{pairwise=FALSE}.
#' If \code{"all"}, the intersection of pairwise marker sets is used, while if \code{"any"}, the union is used.
#'
#' @return If \code{pairwise=TRUE}, a \linkS4class{List} of Lists of character vectors is returned.
#' Each element of the outer list corresponds to cluster X, each element of the inner list corresponds to another cluster Y,
#' and each character vector specifies the marker genes that distinguish X from Y.
#'
#' If \code{pairwise=FALSE}, a List of character vectors is returned.
#' Each character vector contains the marker genes that distinguish X from any/all other clusters. 
#'
#' @author Aaron Lun
#'
#' @details
#' This is a convenience utility that converts the results of pairwise comparisons into a marker list
#' that can be used in downstream functions, e.g., as the marker sets in \pkg{SingleR}.
#' 
#' @seealso
#' \code{\link{pairwiseTTests}} and friends, to obtain \code{de.lists} and \code{pairs}.
#'
#' \code{\link{combineMarkers}}, for another function that consolidates pairwise DE comparisons.
#'
#' @examples
#' data(example.sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(example.sce)), centers=3) 
#' 
#' out <- pairwiseTTests(logcounts(example.sce), 
#'      clusters=paste0("Cluster", kout$cluster))
#' 
#' top <- getTopMarkers(out$statistics, out$pairs)
#' top[[1]]
#' top[[1]][[2]]
#' 
#' top <- getTopMarkers(out$statistics, out$pairs, pairwise=FALSE)
#' top[[1]]
#' @export
#' @importFrom S4Vectors List
#' @importFrom utils head
getTopMarkers <- function(de.lists, pairs, n=10, field="p.value", pairwise=TRUE, combine.type=c("any", "all")) {
    markers <- List()
    all.labels <- sort(unique(c(pairs$first, pairs$second)))
    combine.type <- match.arg(combine.type)

    for (first in all.labels) {
        cur.markers <- List()
        for (second in all.labels) {
            chosen <- which(pairs$first==first & pairs$second==second)
            if (!length(chosen)) {
                cur.markers[[second]] <- character(0)
            } else if (length(chosen)!=1L){ 
                stop(sprintf("multiple entries in 'pairs' for '%s' vs '%s'", first, second))
            } else {
                cur.stats <- de.lists[[chosen]]
                o <- order(cur.stats[[field]])
                cur.markers[[second]] <- rownames(cur.stats)[head(o, n)]
            }
        }

        if (!pairwise) {
            FUN <- if (combine.type=="any") union else intersect
            cur.markers <- Reduce(FUN, cur.markers)
        }
        markers[[first]] <- cur.markers
    }

    markers
}
