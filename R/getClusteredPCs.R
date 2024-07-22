#' Use clusters to choose the number of PCs
#'
#' Cluster cells after using varying number of PCs,
#' and pick the number of PCs using a heuristic based on the number of clusters.
#'
#' @param pcs A numeric matrix of PCs, where rows are cells and columns are dimensions representing successive PCs.
#' @param FUN A clustering function that takes a numeric matrix with rows as cells and
#' returns a vector containing a cluster label for each cell.
#' Defaults to \code{\link{clusterRows}}.
#' @param ... Further arguments to pass to \code{FUN}.
#' Ignored if \code{FUN=NULL}, use \code{BLUSPARAM} instead.
#' @param BLUSPARAM A \linkS4class{BlusterParam} object specifying the clustering to use when \code{FUN=NULL}.
#' @param min.rank Integer scalar specifying the minimum number of PCs to use.
#' @param max.rank Integer scalar specifying the maximum number of PCs to use.
#' @param by Integer scalar specifying what intervals should be tested between \code{min.rank} and \code{max.rank}.
#'
#' @return A \linkS4class{DataFrame} with one row per tested number of PCs.
#' This contains the fields:
#' \describe{
#' \item{\code{n.pcs}:}{Integer scalar specifying the number of PCs used.}
#' \item{\code{n.clusters}:}{Integer scalar specifying the number of clusters identified.}
#' \item{\code{clusters}:}{A \linkS4class{List} containing the cluster identities for this number of PCs.}
#' }
#' The metadata of the DataFrame contains \code{chosen}, 
#' an integer scalar specifying the \dQuote{ideal} number of PCs to use.
#'
#' @details
#' Assume that the data contains multiple subpopulations, each of which is separated from the others on a different axis.
#' For example, each subpopulation could be defined by a unique set of marker genes that drives separation on its own PC.
#' If we had \eqn{x} subpopulations, we would need at least \eqn{x-1} PCs to successfully distinguish all of them.
#' This motivates the choice of the number of PCs provided we know the number of subpopulations in the data.
#'
#' In practice, we do not know the number of subpopulations so we use the number of clusters as a proxy instead.
#' We apply a clustering function \code{FUN} on the first \eqn{d} PCs,
#' and only consider the values of \eqn{d} that yield no more than \eqn{d+1} clusters.
#' If we see more clusters with fewer dimensions, 
#' we consider this to represent overclustering rather than distinct subpopulations,
#' as multiple subpopulations should not be distinguishable on the same axes (based on the assumption above).
#'
#' We choose \eqn{d} that satisfies the constraint above and maximizes the number of clusters.
#' The idea is that more PCs should include more biological signal, allowing \code{FUN} to detect more distinct subpopulations;
#' until the point that the extra signal outweights the added noise at high dimensions,
#' such that resolution decreases and it becomes more difficult for \code{FUN} to distinguish between subpopulations.
#' 
#' Any \code{FUN} can be used that automatically chooses the number of clusters based on the data.
#' The default is a graph-based clustering method using \code{\link{makeSNNGraph}} and \code{\link{cluster_walktrap}},
#' where arguments in \code{...} are passed to the former.
#' Users should not supply \code{FUN} where the number of clusters is fixed in advance, 
#' (e.g., k-means, hierarchical clustering with known \code{k} in \code{\link{cutree}}).
#'
#' The identities of the output clusters are returned at each step for comparison, e.g., using methods like \pkg{clustree}.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' sce <- scater::runPCA(sce)
#' output <- getClusteredPCs(reducedDim(sce))
#' output
#' 
#' metadata(output)$chosen
#' 
#' @seealso
#' \code{\link{runPCA}}, to compute the PCs in the first place.
#'
#' \code{\link{clusterRows}} and \linkS4class{BlusterParam}, for possible choices of \code{BLUSPARAM}.
#'
#' @export
#' @importFrom igraph cluster_walktrap
#' @importFrom S4Vectors DataFrame metadata metadata<- List
#' @importFrom bluster NNGraphParam
getClusteredPCs <- function(pcs, FUN=NULL, ..., BLUSPARAM=NNGraphParam(), min.rank=5, max.rank=ncol(pcs), by=1) {
    if (is.null(FUN)) {
        if (length(list(...))) {
            warning("arguments in '...' are now ignored with 'FUN=NULL'")
        }
        FUN <- function(x, ...) clusterRows(x, BLUSPARAM=BLUSPARAM)
    }

    max.rank <- max(max.rank, min.rank)
    n <- seq(as.integer(min.rank), as.integer(max.rank), by=as.integer(by))
    collected <- vector("list", length(n))
    for (i in seq_along(n)) {
        collected[[i]] <- FUN(pcs[,seq_len(n[i]),drop=FALSE], ...)
    }

    nc <- vapply(collected, FUN=function(x) length(unique(x)), 0L)
    output <- DataFrame(n.pcs=n, n.clusters=nc, clusters=I(List(collected)))

    keep <- which(nc <= n+1)
    if (length(keep)) {
        chosen <- n[keep][which.max(nc[keep])]
    } else {
        chosen <- max.rank # keep the max, that most satisfies the constraint.
    }
    metadata(output)$chosen <- chosen

    output
}
