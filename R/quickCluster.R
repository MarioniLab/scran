setGeneric("quickCluster", function(x, ...) standardGeneric("quickCluster"))

setMethod("quickCluster", "matrix", function(x, min.size=200, subset.row=NULL, ...)  
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach, with modifications by Aaron Lun
# created 1 December 2015
# last modified 5 June 2016
{   
    if (ncol(x) < min.size){
        stop('fewer cells than the mininimum cluster size')
    }
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    distM <- .Call(cxx_compute_cordist, x, subset.row - 1L) # taken into C++ to improve memory efficiency.
    if (is.character(distM)) { stop(distM) }
    distM <- as.dist(distM)
    htree <- hclust(distM, method='ward.D2')
    clusters <- unname(cutreeDynamic(htree, minClusterSize=min.size, distM=as.matrix(distM), verbose=0, ...))

    unassigned <- clusters==0L
    if (any(unassigned)) { 
        warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
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

