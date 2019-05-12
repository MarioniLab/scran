#' @importFrom stats hclust dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom scater librarySizeFactors normalizeCounts
#' @importFrom igraph cluster_walktrap
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom Matrix t
#' @importFrom BiocSingular ExactParam bsdeferred
#' @importClassesFrom Matrix dgCMatrix
.quick_cluster <- function(x, min.size=100, method=c("igraph", "hclust"), use.ranks=NULL,
    d=NULL, subset.row=NULL, min.mean=1, graph.fun=cluster_walktrap,
    BSPARAM=ExactParam(), BPPARAM=SerialParam(), block=NULL, block.BPPARAM=SerialParam(), 
    ...)
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
            MoreArgs=list(min.size=min.size, method=method, use.ranks=use.ranks,
                subset.row=subset.row, min.mean=min.mean, 
                BSPARAM=BSPARAM, BPPARAM=BPPARAM, ...), 
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
    if (is.null(use.ranks)) {
        .Deprecated(msg="Setting 'use.ranks=TRUE' for the old defaults.\nSet 'use.ranks=FALSE' for the new defaults.") 
        use.ranks <- TRUE
    }

    if (use.ranks) {
        if (is.null(d)) {
            d <- 50
        }

        # Only preserving sparsity in the ranks if:
        # - the input type is amenable,
        # - we're using the output for PCA via fast %*%,
        # - we're allowing deferred calculations for PCA.
        deferred <- !is.na(d) && is(x, "dgCMatrix") && bsdeferred(BSPARAM)

        y <- .create_rank_matrix(x, deferred=deferred, subset.row=subset.row, min.mean=min.mean)
    } else {
        sf <- librarySizeFactors(x, subset_row=subset.row)
        y <- normalizeCounts(x, size_factors=sf, return_log=TRUE, subset_row=subset.row)
        if (is.null(d)) {
            fit <- trendVar(y)
            y <- denoisePCA(y, technical=fit$trend, BSPARAM=BSPARAM)
            d <- NA
        } else {
            y <- t(y)
        }
    }

    if (!is.na(d)) {
        svd.out <- .centered_SVD(y, max.rank=d, keep.right=FALSE, BSPARAM=BSPARAM)
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

#' @importFrom Matrix colMeans t
#' @importFrom BiocSingular DeferredMatrix
.create_rank_matrix <- function(x, deferred, ...) {
    if (!deferred) {
        y <- scaledColRanks(x, ..., transposed=TRUE)
    } else {
        y <- scaledColRanks(x, ..., transposed=FALSE, as.sparse=TRUE)
        y <- t(DeferredMatrix(y, center=colMeans(y)))
    }
    y
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

    as.integer(factor(clusters))
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

