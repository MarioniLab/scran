#' Use clusters to choose the number of PCs
#'
#' Cluster cells after using varying number of PCs,
#' and pick the number of PCs using a heuristic based on the number of clusters.
#'
#' @param pcs A numeric matrix of PCs, where rows are cells and columns are dimensions representing successive PCs.
#' @param FUN A clustering function that takes a numeric matrix with rows as cells and
#' returns a vector containing a cluster label for each cell.
#' @param ... Further arguments to pass to \code{FUN}.
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
#' Any \code{FUN} can be used that automatically chooses the number of clusters based on the data.
#' The default is a graph-based clustering method using \code{\link{buildSNNGraph}} and \code{\link{cluster_walktrap}},
#' where arguments in \code{...} are passed to the former.
#' Users should not supply \code{FUN} where the number of clusters is fixed in advance, 
#' (e.g., k-means, hierarchical clustering with known \code{k} in \code{\link{cutree}}).
#' 
#' The idea is that more PCs should include more biological signal, allowing \code{FUN} to detect more distinct clusters.
#' However, this only holds up to the point that the extra signal outweights the added noise at high dimensions.
#' With too many PCs, we would expect to see a decrease in the number of clusters 
#' as it becomes more difficult for \code{FUN} to distinguish between them.
#' The \dQuote{ideal} number of PCs can thus be chosen at the maximum number of clusters.
#'
#' We add another constraint that the number of clusters must be no greater than the number of PCs minus 1.
#' This is because we need at least \code{d} PCs to guarantee that \code{d+1} subpopulations are separated.
#' (For example, in the most extreme case, each subpopulation could be defined by a unique set of marker genes driving its own PC.)
#' The aim is to avoid situations where \code{FUN} generates many clusters at low rank due to the near-absence of noise.
#' This would retain very few PCs and may discard more subtle biological factors.
#'
#' The identities of the output clusters are returned at each step for use in packages like \pkg{clustree}.
#'
#' @author Aaron Lun
#'
#' @examples
#' sce <- scater::mockSCE()
#' sce <- scater::logNormCounts(sce)
#' sce <- scater::runPCA(sce)
#' 
#' output <- getClusteredPCs(reducedDim(sce))
#' output
#' 
#' metadata(output)$chosen
#' 
#' @seealso
#' \code{\link{runPCA}}, to compute the PCs in the first place.
#'
#' \code{\link{buildSNNGraph}}, for arguments to use in with default \code{FUN}.
#'
#' @export
#' @importFrom igraph cluster_walktrap
#' @importFrom S4Vectors DataFrame metadata metadata<- List
getClusteredPCs <- function(pcs, FUN=NULL, ..., min.rank=5, max.rank=ncol(pcs), by=1) {
    if (is.null(FUN)) {
        FUN <- function(x, ...) {
            g <- buildSNNGraph(x, ..., transposed=TRUE)
            cluster_walktrap(g)$membership
        }
    }

    max.rank <- max(max.rank, min.rank)
    n <- seq(as.integer(min.rank), as.integer(max.rank), by=as.integer(by))
    collected <- vector("list", length(n))
    for (i in seq_along(n)) {
        collected[[i]] <- FUN(pcs[,seq_len(n[i]),drop=FALSE], ...)
    }

    nc <- vapply(collected, FUN=function(x) length(unique(x)), 0L)
    output <- DataFrame(n.pcs=n, n.clusters=nc, clusters=I(List(collected)))

    keep <- which(nc < n-1)
    if (length(keep)) {
        chosen <- n[keep][which.max(nc[keep])]
    } else {
        chosen <- max.rank # keep the max, that most satisfies the constraint.
    }
    metadata(output)$chosen <- chosen

    output
}
