.quick_cluster <- function(x, min.size=200, max.size=NULL, subset.row=NULL, get.ranks=FALSE, 
                           method=c("hclust", "igraph"), pc.approx=TRUE, ...)  
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach, with modifications by Aaron Lun
# created 1 December 2015
# last modified 3 October 2017
{   
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

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
    } else if (!is.null(max.size) && max.size < min.size*2) { 
        stop("maximum cluster size must be at least twice the minimum size") 
        # otherwise, split clusters aren't guaranteed to be > min.size. See below.
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

    clusters <- .limit_cluster_size(clusters, max.size)
    clusters <- factor(clusters)
    return(clusters)
}

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
            next.m <- modularity(g, next.clust)
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

.limit_cluster_size <- function(clusters, max.size) {
    if (is.null(max.size)) { return(clusters) }
    
    new.clusters <- clusters
    counter <- 1L
    for (id in unique(clusters)) {
        current <- id==clusters
        ncells <- sum(current)
        
        if (ncells <= max.size) {
            new.clusters[current] <- counter
            counter <- counter+1L
            next
        }
       
        # Size of output clusters is max.size * N / ceil(N), where N = ncells/max.size.
        # This is minimal at the smallest N > 1, where output clusters are at least max.size/2. 
        # Thus, we need max.size/2 >= min.size to guarantee that the output clusters are >= min.size.
        mult <- ceiling(ncells/max.size)
        realloc <- rep(seq_len(mult) - 1L + counter, length.out=ncells)
        new.clusters[current] <- realloc
        counter <- counter + mult
    }
    return(new.clusters)
}

setGeneric("quickCluster", function(x, ...) standardGeneric("quickCluster"))

setMethod("quickCluster", "ANY", .quick_cluster)

setMethod("quickCluster", "SingleCellExperiment", 
          function(x, subset.row=NULL, ..., assay.type="counts", get.spikes=FALSE) { 

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)          
    .quick_cluster(assay(x, i=assay.type), subset.row=subset.row, ...)
})

