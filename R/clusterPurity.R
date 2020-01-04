#' Evaluate cluster purity
#'
#' Use a hypersphere-based approach to compute the \dQuote{purity} of each cluster based on the number of contaminating cells in its region of the coordinate space.
#'
#' @param x For the ANY method, a numeric matrix-like object containing expression values for genes (rows) and cells (columns).
#' If \code{transposed=TRUE}, cells should be in rows and reduced dimensions should be in the columns.
#'
#' For the \linkS4class{SingleCellExperiment} method, a SingleCellExperiment object containing an expression matrix.
#' If \code{use.dimred} is supplied, it should contain a reduced dimension result in its \code{\link{reducedDims}}.
#' @param clusters Factor specifying the cluster identity for each cell.
#' @param k Integer scalar specifying the number of nearest neighbors to use to determine the radius of the hyperspheres.
#' @param transposed Logical scalar specifying whether \code{x} contains cells in the rows.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to use for purity calculations.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' This should be an algorithm supported by \code{\link{findNeighbors}}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes.
#'
#' @return
#' A numeric vector of purity values for each cell in \code{x}.
#' 
#' @details
#' The purity of a cluster is quantified by creating a hypersphere around each cell in the cluster
#' and computing the proportion of cells in that hypersphere from the same cluster.
#' If all cells in a cluster have proportions close to 1, this indicate that the cluster is highly pure,
#' i.e., there are no cells from other clusters in its region of the coordinate space.
#'
#' A more common scenario is that a fraction of cells will have low purity values;
#' this corresponds to cells lying at the boundaries of two adjacent clusters,
#' whereas cells close to the cluster center will have relatively high purity.
#' A high degree of over-clustering will manifest as a majority of cells with purities close to zero.
#' 
#' The choice of \code{k} is used only to determine an appropriate value for the hypersphere radius.
#' We use hyperspheres as this is robust to changes in density throughout the coordinate space,
#' in contrast to computing purity based on the proportion of k-nearest neighbors in the same cluster.
#' For example, the latter will fail most obviously when the size of the cluster is less than \code{k}.
#'
#' Technically, purity values are computed by adding \code{pseudo.count} to the number of cells in the same cluster
#' and dividing by the total number of cells in the hypersphere plus \code{2*pseudo.count}.
#' This shrinks all values towards 0.5 and will yield a more conservative estimate of purity when few cells are available.
#'
#' @author Aaron Lun
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' g <- buildSNNGraph(sce)
#' clusters <- igraph::cluster_walktrap(g)$membership
#' out <- clusterPurity(sce, clusters)
#' boxplot(split(out, clusters))
#'
#' # Mocking up a stronger example:
#' ngenes <- 1000
#' centers <- matrix(rnorm(ngenes*3), ncol=3)
#' clusters <- sample(1:3, 500, replace=TRUE)
#'
#' y <- centers[,clusters]
#' y <- y + rnorm(length(y))
#' 
#' out2 <- clusterPurity(y, clusters)
#' boxplot(split(out2, clusters))
#'
#' @name clusterPurity
NULL

#' @importFrom BiocNeighbors KmknnParam buildIndex findKNN findNeighbors
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scater .bpNotSharedOrUp
#' @importFrom Matrix t
#' @importFrom stats median
.cluster_purity <- function(x, clusters, k=50, transposed=FALSE, subset.row=NULL,
    pseudo.count=1, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    }

    x <- as.matrix(x)
    idx <- buildIndex(x, BNPARAM=BNPARAM)

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    dist <- median(findKNN(BNINDEX=idx, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, last=1, get.index=FALSE)$distance)
    nout <- findNeighbors(BNINDEX=idx, threshold=dist, BNPARAM=BNPARAM, get.distance=FALSE)$index

    # Counting the number of entities in the same cluster as me.
    # This could be made much more efficient.
    sameness <- clusters[unlist(nout)]==rep(clusters, lengths(nout))
    samelist <- relist(sameness, nout)
    nsame <- vapply(samelist, sum, 0L)

    (nsame + pseudo.count)/(lengths(samelist) + 2*pseudo.count)
}

#' @export
#' @rdname clusterPurity
setGeneric("clusterPurity", function(x, ...) standardGeneric("clusterPurity"))

#' @export
#' @rdname clusterPurity
setMethod("clusterPurity", "ANY", .cluster_purity)

#' @export
#' @rdname clusterPurity
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
setMethod("clusterPurity", "SingleCellExperiment", function(x, ..., assay.type="logcounts", use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        transposed <- TRUE
        x <- reducedDim(x, use.dimred)
    } else {
        x <- assay(x, assay.type)
        transposed <- FALSE
    }
    .cluster_purity(x, ..., transposed=transposed)
})
