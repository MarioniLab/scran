#' @importFrom stats hclust dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom scater calcAverage
#' @importFrom igraph cluster_fast_greedy
#' @importFrom BiocParallel SerialParam bpmapply
.quick_cluster <- function(x, min.size=200, max.size=NULL, method=c("hclust", "igraph"),
                           pc.approx=TRUE, get.ranks=FALSE, subset.row=NULL, min.mean=1, 
                           block=NULL, block.BPPARAM=SerialParam(), ...)
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach, with modifications by Aaron Lun
# created 1 December 2015
{
    if (!is.null(block) && length(unique(block))>1L) {
        # Splitting into parallel processes across blocks.
        # We create submatrices here to avoid memory allocation within each core.
        by.block <- split(seq_along(block), block)
        x.by.block <- .split_matrix_by_workers(x, by.block, byrow=FALSE)
        collected <- bpmapply(FUN=.quick_cluster, x.by.block, 
            MoreArgs=list(min.size=min.size, max.size=max.size, method=method, 
                pc.approx=pc.approx, get.ranks=get.ranks, subset.row=subset.row, 
                min.mean=min.mean, ...), 
            BPPARAM=block.BPPARAM)

        # Merging the results across different blocks.
        reordering <- order(unlist(by.block, use.names=FALSE))
        if (get.ranks) {
            return(do.call(cbind, collected)[,reordering,drop=FALSE])
        } else {
            last <- 0L
            for (b in seq_along(collected)) {
                to.add <- nlevels(collected[[b]])
                collected[[b]] <- as.integer(collected[[b]]) + last
                last <- last + to.add
            }
            return(factor(unlist(collected, use.names=FALSE)[reordering]))
        }
    }

    rkout <- scaledColRanks(x, subset.row=subset.row, min.mean=min.mean, transposed=!get.ranks)
    if (get.ranks) {
        .Deprecated(msg="'get.ranks=TRUE' is deprecated.\nUse 'scaledColRanks' instead.")
        return(rkout)
    }

    # Checking size specifications.
    if (ncol(x) < min.size){
        stop('fewer cells than the minimum cluster size')
    } 

    method <- match.arg(method)
    if (method=="igraph") { 
        g <- buildSNNGraph(rkout, pc.approx=pc.approx, transposed=TRUE, ...)
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

