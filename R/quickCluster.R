setGeneric("quickCluster", function(x, ...) standardGeneric("quickCluster"))

setMethod("quickCluster", "matrix", function(x, min.size=200, subset.row=NULL, get.ranks=FALSE, method=c("hclust", "igraph"), ...)  
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach, with modifications by Aaron Lun
# created 1 December 2015
# last modified 4 April 2017
{   
    if (ncol(x) < min.size){
        stop('fewer cells than the minimum cluster size')
    }
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    get.ranks <- as.logical(get.ranks)

    method <- match.arg(method)
    rkout <- .Call(cxx_compute_cordist, x, subset.row - 1L, get.ranks | method=="igraph") # taken into C++ to improve memory efficiency.
    if (is.character(rkout)) { stop(rkout) }
    if (get.ranks) {
        return(rkout)
    }

    if (method=="igraph") { 
        g <- buildSNNGraph(rkout, ...)
        out <- cluster_fast_greedy(g)
        clusters <- out$membership

    } else {
        distM <- as.dist(rkout)
        htree <- hclust(distM, method='ward.D2')
        clusters <- unname(cutreeDynamic(htree, minClusterSize=min.size, distM=rkout, verbose=0, ...))
        
        unassigned <- clusters==0L
        if (any(unassigned)) { 
            warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
        }
    }
    clusters <- factor(clusters)
    return(clusters)
})

setMethod("quickCluster", "SCESet", function(x, subset.row=NULL, ..., assay="counts", get.spikes=FALSE) { 
    if (is.null(subset.row)) {
        subset.row <- .spikeSubset(x, get.spikes)
    }
    quickCluster(assayDataElement(x, assay), subset.row=subset.row, ...)
})

