#' Quick clustering of cells
#'
#' Cluster similar cells based on their expression profiles, using either log-expression values or ranks.
#' 
#' @param x A numeric count matrix where rows are genes and columns are cells.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param min.size An integer scalar specifying the minimum size of each cluster.
#' @param method String specifying the clustering method to use.
#' \code{"hclust"} uses hierarchical clustering while \code{"igraph"} uses graph-based clustering.
#' @param use.ranks A logical scalar indicating whether clustering should be performed on the rank matrix, 
#' i.e., based on Spearman's rank correlation.
#' @param d An integer scalar specifying the number of principal components to retain.
#' If \code{d=NULL} and \code{use.ranks=TRUE}, this defaults to 50.
#' If \code{d=NULL} and \code{use.rank=FALSE}, the number of PCs is chosen by \code{\link{denoisePCA}}.
#' If \code{d=NA}, no dimensionality reduction is performed and the gene expression values (or their rank equivalents) are directly used in clustering.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param min.mean A numeric scalar specifying the filter to be applied on the average count for each filter prior to computing ranks.
#' Only used when \code{use.ranks=TRUE}, see \code{?\link{scaledColRanks}} for details.
#' @param graph.fun A function specifying the community detection algorithm to use on the nearest neighbor graph when \code{method="igraph"}.
#' Usually obtained from the \pkg{igraph} package.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, if \code{d} is not \code{NA}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object to use for parallel processing within each block.
#' @param block A factor of length equal to \code{ncol(x)} specifying whether clustering should be performed within pre-specified blocks.
#' By default, all columns in \code{x} are treated as a single block.
#' @param block.BPPARAM A \linkS4class{BiocParallelParam} object specifying whether and how parallelization should be performed across blocks, 
#' if \code{block} is non-\code{NULL} and has more than one level.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY method, additional arguments to be passed to \code{\link{NNGraphParam}} for \code{method="igraph"},
#' or to be included in \code{cut.params} argument for \code{\link{HclustParam}} when \code{method="hclust"}.
#' 
#' For the \linkS4class{SummarizedExperiment} method, additional arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use. 
#' 
#' @details
#' This function provides a convenient wrapper to quickly define clusters of a minimum size \code{min.size}.
#' Its intended use is to generate \dQuote{quick and dirty} clusters for use in \code{\link{computeSumFactors}}.
#' Two clustering strategies are available:
#' \itemize{
#' \item If \code{method="hclust"}, a distance matrix is constructed;
#' hierarchical clustering is performed using Ward's criterion;
#' and \code{\link[dynamicTreeCut]{cutreeDynamic}} is used to define clusters of cells.
#' \item If \code{method="igraph"}, a shared nearest neighbor graph is constructed using the \code{\link{buildSNNGraph}} function.
#' This is used to define clusters based on highly connected communities in the graph, using the \code{graph.fun} function.
#' }
#' 
#' By default, \code{quickCluster} will apply these clustering algorithms on the principal component (PC) scores generated from the log-expression values.
#' These are obtained by running \code{\link{denoisePCA}} on HVGs detected using the trend fitted to endogenous genes with \code{\link{modelGeneVar}}.
#' If \code{d} is specified, the PCA is directly performed on the entire \code{x} and the specified number of PCs is retained.
#'
#' % We use modelGeneVar because we don't know if there's spike-ins or if it's UMI data or not.
#' % This could behave weirdly with denoisePCA, but it given that this function is intended for quick-and-dirty clustering,
#' % some inappropriate loss of signal is not criticial.
#'
#' It is also possible to use the clusters from this function for actual biological interpretation.
#' In such cases, users should set \code{min.size=0} to avoid aggregation of small clusters.
#' However, it is often better to call the relevant functions (\code{\link{modelGeneVar}}, \code{\link{denoisePCA}} and \code{\link{buildSNNGraph}}) manually as this provides more opportunities for diagnostics when the meaning of the clusters is important.
#'
#' @section Clustering within blocks:
#' We can break up the dataset by specifying \code{block} to cluster cells, usually within each batch or run.
#' This generates clusters within each level of \code{block}, which is entirely adequate for applications like \code{\link{computeSumFactors}} where the aim of clustering is to separate dissimilar cells rather than group together all similar cells.
#' Blocking reduces computational work considerably by allowing each level to be processed independently, without compromising performance provided that there are enough cells within each batch.
#' 
#' Indeed, for applications like \code{\link{computeSumFactors}}, we can use \code{block} even in the absence of any known batch structure.
#' Specifically, we can set it to an arbitrary factor such as \code{block=cut(seq_len(ncol(x)), 10)} to split the cells into ten batches of roughly equal size.
#' This aims to improve speed, especially when combined with \code{block.PARAM} to parallelize processing of the independent levels.
#'
#' @section Using ranks:
#' If \code{use.ranks=TRUE}, clustering is instead performed on PC scores obtained from scaled and centred ranks generated by \code{\link{scaledColRanks}}.
#' This effectively means that clustering uses distances based on the Spearman's rank correlation between two cells.
#' In addition, if \code{x} is a \linkS4class{dgCMatrix} and \code{BSPARAM} has \code{deferred=TRUE}, 
#' ranks will be computed without loss of sparsity to improve speed and memory efficiency during PCA.
#'
#' When \code{use.ranks=TRUE}, the function will filter out genes with average counts (as defined by \code{\link{calculateAverage}}) below \code{min.mean} prior to computing ranks.
#' This removes low-abundance genes with many tied ranks, especially due to zeros, which may reduce the precision of the clustering.
#' We suggest setting \code{min.mean} to 1 for read count data and 0.1 for UMI data - the function will automatically try to determine this from the data if \code{min.mean=NULL}.
#' 
#' Setting \code{use.ranks=TRUE} is invariant to scaling normalization and avoids circularity between normalization and clustering, e.g., in \code{\link{computeSumFactors}}.
#' However, the default is to use the log-expression values with \code{use.ranks=FALSE}, as this yields finer and more precise clusters.
#' 
#' @section Enforcing cluster sizes:
#' With \code{method="hclust"}, \code{\link[dynamicTreeCut]{cutreeDynamic}} is used to ensure that all clusters contain a minimum number of cells.
#' However, some cells may not be assigned to any cluster and are assigned identities of \code{"0"} in the output vector.
#' In most cases, this is because those cells belong in a separate cluster with fewer than \code{min.size} cells.
#' The function will not be able to call this as a cluster as the minimum threshold on the number of cells has not been passed.
#' Users are advised to check that the unassigned cells do indeed form their own cluster.
#' Otherwise, it may be necessary to use a different clustering algorithm.
#' 
#' When using \code{method="igraph"}, clusters are first identified using the specified \code{graph.fun}.
#' If the smallest cluster contains fewer cells than \code{min.size}, it is merged with the closest neighbouring cluster.
#' In particular, the function will attempt to merge the smallest cluster with each other cluster.
#' The merge that maximizes the modularity score is selected, and a new merged cluster is formed.
#' This process is repeated until all (merged) clusters are larger than \code{min.size}.
#' 
#' @return
#' A character vector of cluster identities for each cell in \code{x}.
#' 
#' @author
#' Aaron Lun and Karsten Bach
#' 
#' @seealso
#' \code{\link{computeSumFactors}}, where the clustering results can be used as \code{clusters=}.
#' 
#' \code{\link{buildSNNGraph}}, for additional arguments to customize the clustering when \code{method="igraph"}.
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}}, for additional arguments to customize the clustering when \code{method="hclust"}.
#' 
#' \code{\link{scaledColRanks}}, to get the rank matrix that was used with \code{use.rank=TRUE}.
#'
#' \code{\link{quickSubCluster}}, for a related function that uses a similar approach for subclustering.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' 
#' # Basic application (lowering min.size for this small demo):
#' clusters <- quickCluster(sce, min.size=50)
#' table(clusters)
#' 
#' # Operating on ranked expression values:
#' clusters2 <- quickCluster(sce, min.size=50, use.ranks=TRUE)
#' table(clusters2)
#' 
#' # Using hierarchical clustering:
#' clusters <- quickCluster(sce, min.size=50, method="hclust")
#' table(clusters)
#' @references
#' van Dongen S and Enright AJ (2012).
#' Metric distances derived from cosine similarity and Pearson and Spearman correlations.
#' \emph{arXiv} 1208.3145
#' 
#' Lun ATL, Bach K and Marioni JC (2016).
#' Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
#' \emph{Genome Biol.} 17:75
#'
#' @name quickCluster
NULL

#' @importFrom stats hclust dist
#' @importFrom scuttle librarySizeFactors normalizeCounts .guessMinMean
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom Matrix t
#' @importFrom BiocSingular bsparam bsdeferred
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom bluster clusterRows NNGraphParam HclustParam
.quick_cluster <- function(x, min.size=100, method=c("igraph", "hclust"), use.ranks=FALSE,
    d=NULL, subset.row=NULL, min.mean=NULL, graph.fun="walktrap",
    BSPARAM=bsparam(), BPPARAM=SerialParam(), block=NULL, block.BPPARAM=SerialParam(), 
    ...)
{
    if (!is.null(block) && length(unique(block))>1L) {
        # Splitting into parallel processes across blocks.
        # We create submatrices here to avoid memory allocation within each core.
        by.block <- split(seq_along(block), block)
        x.by.block <- lapply(by.block, function(i) x[,i,drop=FALSE])

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
    if (use.ranks) {
        if (is.null(d)) {
            d <- 50
        }

        # Only preserving sparsity in the ranks if:
        # - the input type is amenable,
        # - we're using the output for PCA via fast %*%,
        # - we're allowing deferred calculations for PCA.
        deferred <- !is.na(d) && is(x, "dgCMatrix") && bsdeferred(BSPARAM)

        min.mean <- .guessMinMean(x, min.mean=min.mean, BPPARAM=BPPARAM)
        y <- .create_rank_matrix(x, deferred=deferred, subset.row=subset.row, min.mean=min.mean, BPPARAM=BPPARAM)
    } else {
        sf <- librarySizeFactors(x, subset_row=subset.row)
        y <- normalizeCounts(x, size_factors=sf, subset_row=subset.row)
        if (is.null(d)) {
            fit <- modelGeneVar(y)
            chosen <- getTopHVGs(fit, n=500, prop=0.1) # At least 500 genes, or 10% of genes; whichever is larger.
            y <- getDenoisedPCs(y, technical=fit, subset.row=chosen, BSPARAM=BSPARAM)$components
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
        out <- clusterRows(y, NNGraphParam(..., cluster.fun=graph.fun), full=TRUE)
        clusters <- .merge_closest_graph(out$objects$graph, as.integer(out$clusters), min.size=min.size)
        clusters <- factor(clusters)

    } else {
        clusters <- clusterRows(y, HclustParam(method="ward.D2", cut.dynamic=TRUE, cut.params=list(minClusterSize=min.size, ...)))
        unassigned <- clusters=="0"
        if (any(unassigned)) { 
            warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
        }
    }

    clusters
}

############################
# Internal functions.
############################

#' @importFrom Matrix colMeans t
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom BiocParallel SerialParam
.create_rank_matrix <- function(x, deferred, ..., BPPARAM=SerialParam()) {
    if (!deferred) {
        y <- scaledColRanks(x, ..., transposed=TRUE)
    } else {
        old <- getAutoBPPARAM()
        setAutoBPPARAM(BPPARAM)
        on.exit(setAutoBPPARAM(old))

        y <- scaledColRanks(x, ..., transposed=FALSE, as.sparse=TRUE, BPPARAM=BPPARAM)
        y <- t(ScaledMatrix::ScaledMatrix(y, center=colMeans(y)))
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
#' @rdname quickCluster
setGeneric("quickCluster", function(x, ...) standardGeneric("quickCluster"))

#' @export
#' @rdname quickCluster
setMethod("quickCluster", "ANY", .quick_cluster)

#' @export
#' @rdname quickCluster
#' @importFrom SummarizedExperiment assay
setMethod("quickCluster", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .quick_cluster(assay(x, i=assay.type), ...)
})

