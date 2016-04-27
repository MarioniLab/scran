setGeneric("quickCluster", function(x, ...) standardGeneric("quickCluster"))

setMethod("quickCluster", "matrix", function(x, min.size=200, ...)  
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach, with modifications by Aaron Lun
# created 1 December 2015
# last modified 17 February 2016
{   
    if (ncol(x) < min.size){
        stop('fewer cells than the mininimum cluster size')
    }

    distM <- as.dist( 1 - cor(x, method='spearman'))
    htree <- hclust(distM, method='ward.D2')
    clusters <- unname(cutreeDynamic(htree, minClusterSize=min.size, distM=as.matrix(distM), verbose=0, ...))

    unassigned <- clusters==0L
    if (any(unassigned)) { 
        warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
    }
    clusters <- factor(clusters)
    return(clusters)
})

setMethod("quickCluster", "SCESet", function(x, ..., assay="counts", get.spikes=FALSE) { 
    quickCluster(.getUsedMatrix(x, assay, get.spikes), ...) 
})

