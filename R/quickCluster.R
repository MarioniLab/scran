#' @importFrom stats hclust dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom scater calcAverage librarySizeFactors normalizeCounts
#' @importFrom igraph cluster_walktrap
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom BiocGenerics t
.quick_cluster <- function(x, min.size=100, method=c("igraph", "hclust"), use.ranks=FALSE,
    pc.approx=FALSE, d=NULL, subset.row=NULL, min.mean=1, graph.fun=cluster_walktrap,
    block=NULL, block.BPPARAM=SerialParam(), ...)
# Generates a factor specifying the cluster to which each cell is assigned.
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
            MoreArgs=list(min.size=min.size, method=method, 
                pc.approx=pc.approx, subset.row=subset.row, 
                min.mean=min.mean, ...), 
            BPPARAM=block.BPPARAM)

        # Merging the results across different blocks.
        reordering <- order(unlist(by.block, use.names=FALSE))
        last <- 0L
        for (b in seq_along(collected)) {
            to.add <- nlevels(collected[[b]])
            collected[[b]] <- as.integer(collected[[b]]) + last
            last <- last + to.add
        }
        return(factor(unlist(collected, use.names=FALSE)[reordering]))
    }

    # Obtaining some values to use for clustering.
    if (use.ranks) {
        y <- scaledColRanks(x, subset.row=subset.row, min.mean=min.mean, transposed=TRUE)
        if (is.null(d)) {
            d <- 50
        }
    } else {
        sf <- librarySizeFactors(x, subset_row=subset.row)
        y <- normalizeCounts(x, size_factors=sf, return_log=TRUE, subset_row=subset.row)
        if (is.null(d)) {
            fit <- trendVar(y)
            y <- denoisePCA(y, technical=fit$trend, approximate=pc.approx)
            d <- NA
        } else {
            y <- t(y)
        }
    }

    if (!is.na(d)) {
        svd.out <- .centered_SVD(y, max.rank=d, approximate=pc.approx, keep.right=FALSE)
        y <- .svd_to_pca(svd.out, d, named=FALSE)
    }

    # Checking size specifications.
    if (ncol(x) < min.size){
        stop('fewer cells than the minimum cluster size')
    } 

    method <- match.arg(method)
    if (method=="igraph") { 
        g <- buildSNNGraph(y, d=NA, transposed=TRUE, ...)
        out <- graph.fun(g)
        clusters <- out$membership
        clusters <- .merge_closest_graph(g, clusters, min.size=min.size)

    } else {
        distM <- dist(as.matrix(y)) # Coercing to matrix, if it isn't already.
        htree <- hclust(distM, method='ward.D2')
        clusters <- unname(cutreeDynamic(htree, minClusterSize=min.size, distM=as.matrix(distM), verbose=0, ...))
        
        unassigned <- clusters==0L
        if (any(unassigned)) { 
            warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
        }
    }

    factor(clusters)
}

############################
# Internal functions.
############################

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

############################
# S4 method definitions 
############################

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

