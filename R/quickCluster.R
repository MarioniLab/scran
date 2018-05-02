#' @importFrom stats hclust dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom scater calcAverage
#' @importFrom igraph cluster_fast_greedy
.quick_cluster <- function(x, min.size=200, max.size=NULL, method=c("hclust", "igraph"),
                           pc.approx=TRUE, get.ranks=FALSE, subset.row=NULL, min.mean=1, ...)
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach, with modifications by Aaron Lun
# created 1 December 2015
{   
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    if (!is.null(min.mean) && all(dim(x)>0L)) {
        further.subset <- calcAverage(x, subset_row=subset.row) >= min.mean
        subset.row <- subset.row[further.subset]
    }

    # Obtaining scaled/centred ranks to compute cosine distances.
    # Using this instead of colRanks to support dgCMatrix, HDF5Array objects.
    get.ranks <- as.logical(get.ranks)
    method <- match.arg(method)
    rkout <- .Call(cxx_get_scaled_ranks, x, subset.row-1L, !get.ranks && method!="igraph")
    if (get.ranks) {
        return(rkout)
    }

    # Checking size specifications.
    if (ncol(x) < min.size){
        stop('fewer cells than the minimum cluster size')
    } 

    if (method=="igraph") { 
        g <- buildSNNGraph(rkout, pc.approx=pc.approx, ...)
        out <- cluster_fast_greedy(g)
        clusters <- out$membership
        clusters <- .merge_closest_graph(g, clusters, min.size=min.size)

    } else {
        distM <- dist(as.matrix(rkout)) # Coercing to matrix, if it isn't already.
        htree <- hclust(distM, method='ward.D2')
        clusters <- unname(cutreeDynamic(htree, minClusterSize=min.size, distM=as.matrix(distM), verbose=0, ...))
        
        unassigned <- clusters==0L
        if (any(unassigned)) { 
            warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
        }
    }

    if (!is.null(max.size)) {
        .Deprecated("'max.size' is deprecated, use 'max.cluster.size=' in 'computeSumFactors' instead")
        clusters <- .limit_cluster_size(clusters, max.size)
    }
    factor(clusters)
}

#' @importFrom igraph modularity E
.merge_closest_graph <- function(g, clusters, min.size) {
    repeat {
        all.sizes <- table(clusters)
        if (all(all.sizes >= min.size)) { break }

        clust.names <- as.integer(names(all.sizes))
        if (length(clust.names)==2L) { 
            clusters[] <- 1L
            break
        }
        
        # Picking the smallest cluster and picking the merge with greatest modularity.
        failed <- clust.names[which.min(all.sizes)]
        to.merge <- clusters==failed
        max.m <- 0
        max.clust <- clusters
        
        for (other in clust.names) {
            if (other==failed) { next }
            next.clust <- clusters
            next.clust[to.merge] <- other
            next.m <- modularity(g, next.clust, weights=E(g)$weight)
            if (max.m < next.m) {
                max.m <- next.m 
                max.clust <- next.clust
            }
        }

        clusters <- max.clust
    }

    clusters <- as.integer(factor(clusters))
    return(clusters)
}

#' @export
setGeneric("quickCluster", function(x, ...) standardGeneric("quickCluster"))

#' @export
setMethod("quickCluster", "ANY", .quick_cluster)

#' @importFrom SummarizedExperiment assay
#' @export
setMethod("quickCluster", "SingleCellExperiment", 
          function(x, subset.row=NULL, ..., assay.type="counts", get.spikes=FALSE) { 

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)          
    .quick_cluster(assay(x, i=assay.type), subset.row=subset.row, ...)
})

