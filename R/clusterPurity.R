#' Evaluate cluster purity
#'
#' Use a hypersphere-based approach to compute the \dQuote{purity} of each cluster based on the number of contaminating cells in its region of the coordinate space.
#'
#' @param x A numeric matrix-like object containing expression values for genes (rows) and cells (columns).
#' If \code{transposed=TRUE}, cells should be in rows and reduced dimensions should be in the columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#'
#' For the SingleCellExperiment method, if \code{use.dimred} is supplied, 
#' \code{x} contain an appropriate reduced dimension result in its \code{\link{reducedDims}}.
#' @param clusters Vector of length equal to \code{ncol(x)}, specifying the cluster identity for each cell.
#' @param k Integer scalar specifying the number of nearest neighbors to use to determine the radius of the hyperspheres.
#' @param transposed Logical scalar specifying whether \code{x} contains cells in the rows.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' This should be an algorithm supported by \code{\link{findNeighbors}}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes.
#' @param weighted A logical scalar indicating whether to weight each cell in inverse proportion to the size of its cluster.
#' Alternatively, a numeric vector of length equal to \code{clusters} containing the weight to use for each cell.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, arguments to pass to the SummarizedExperiment method.
#' @param assay.type A string specifying which assay values to use.
#' @param use.dimred A string specifying whether existing values in \code{reducedDims(x)} should be used.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per cell in \code{x} and the columns:
#' \itemize{
#' \item \code{purity}, a numeric field containing the purity value for the current cell.
#' \item \code{maximum}, the cluster with the highest proportion of cells neighboring the current cell.
#' }
#' Row names are defined as the column names of \code{x}.
#' 
#' @details
#' The purity of a cluster is quantified by creating a hypersphere around each cell in the cluster
#' and computing the proportion of cells in that hypersphere from the same cluster.
#' If all cells in a cluster have proportions close to 1, this indicates that the cluster is highly pure,
#' i.e., there are few cells from other clusters in its region of the coordinate space.
#' The distribution of purities for each cluster can be used as a measure of separation from other clusters.
#'
#' In most cases, the majority of cells of a cluster will have high purities, corresponding to cells close to the cluster center.
#' A fraction of cells will have low values as these lie at the boundaries of two adjacent clusters.
#' A high degree of over-clustering will manifest as a majority of cells with purities close to zero.
#' The \code{maximum} field in the output can be used to determine the identity of the cluster 
#' with the greatest presence in a cell's neighborhood, usually an adjacent cluster for cells lying on the boundary.
#' 
#' The choice of \code{k} is used only to determine an appropriate value for the hypersphere radius.
#' We use hyperspheres as this is robust to changes in density throughout the coordinate space,
#' in contrast to computing purity based on the proportion of k-nearest neighbors in the same cluster.
#' For example, the latter will fail most obviously when the size of the cluster is less than \code{k}.
#'
#' @section Weighting by frequency:
#' By default, purity values are computed after weighting each cell by the reciprocal of the number of cells in the same cluster.
#' Otherwise, clusters with more cells will have higher purities as any contamination is offset by the bulk of cells,
#' which would compromise comparisons of purities between clusters.
#' One can interpret the weighted purities as the expected value after downsampling all clusters to the same size.
#'
#' Advanced users can achieve greater control by manually supplying a numeric vector of weights to \code{weighted}.
#' For example, we may wish to check the purity of batches after batch correction.
#' In this application, \code{clusters} should be set to the \emph{blocking factor} (not the cluster identities!)
#' and \code{weighted} should be set to 1 over the frequency of each combination of cell type and batch.
#' This accounts for differences in cell type composition between batches when computing purities.
#' 
#' If \code{weighted=FALSE}, no weighting is performed.
#'
#' @author Aaron Lun
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' g <- buildSNNGraph(sce)
#' clusters <- igraph::cluster_walktrap(g)$membership
#' out <- clusterPurity(sce, clusters)
#' boxplot(split(out$purity, clusters))
#'
#' # Mocking up a stronger example:
#' ngenes <- 1000
#' centers <- matrix(rnorm(ngenes*3), ncol=3)
#' clusters <- sample(1:3, ncol(sce), replace=TRUE)
#'
#' y <- centers[,clusters]
#' y <- y + rnorm(length(y))
#' 
#' out2 <- clusterPurity(y, clusters)
#' boxplot(split(out2$purity, clusters))
#'
#' @name clusterPurity
NULL

#' @importFrom BiocNeighbors KmknnParam buildIndex findKNN findNeighbors
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom Matrix t
#' @importFrom stats median
#' @importFrom S4Vectors List DataFrame
#' @importClassesFrom IRanges IntegerList LogicalList
#' @importMethodsFrom BiocGenerics relist
.cluster_purity <- function(x, clusters, k=50, transposed=FALSE, weighted=TRUE, 
    subset.row=NULL, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
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

    nout <- List(nout)
    clust.ids <- clusters[unlist(nout)]

    # Constructing weights.
    if (isFALSE(weighted)) {
        w <- rep(1, length(clust.ids))
    } else if (isTRUE(weighted)) {
        w <- 1/table(clusters)
        w <- as.numeric(w[as.character(clust.ids)])
    } else {
        w <- as.numeric(weighted)[unlist(nout)]
    }

    weights <- relist(w, nout)
    total <- sum(weights)

    # Computing the proportions for each cluster.
    uclust <- sort(unique(clusters))
    targets <- matrix(0, length(nout), length(uclust))

    for (i in seq_along(uclust)) {
        current <- clust.ids==uclust[i]
        targets[,i] <- sum(weights[relist(current, nout)])
    }

    same.col <- match(clusters, uclust)

    DataFrame( 
        purity=targets[cbind(seq_along(same.col), same.col)]/total,
        maximum=uclust[max.col(targets, ties.method="first")],
        row.names=rownames(x) # remember, we transposed this earlier.
    )
}

#' @export
#' @rdname clusterPurity
setGeneric("clusterPurity", function(x, ...) standardGeneric("clusterPurity"))

#' @export
#' @rdname clusterPurity
setMethod("clusterPurity", "ANY", .cluster_purity)

#' @export
#' @rdname clusterPurity
#' @importFrom SummarizedExperiment assay
setMethod("clusterPurity", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .cluster_purity(assay(x, assay.type), ..., transposed=FALSE)
})

#' @export
#' @rdname clusterPurity
#' @importFrom SingleCellExperiment reducedDim colLabels
#' @importFrom SummarizedExperiment assay
setMethod("clusterPurity", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"),
    ..., assay.type="logcounts", use.dimred=NULL)
{
    force(clusters)
    if (!is.null(use.dimred)) {
        transposed <- TRUE
        x <- reducedDim(x, use.dimred)
    } else {
        x <- assay(x, assay.type)
        transposed <- FALSE
    }
    .cluster_purity(x, clusters=clusters, ..., transposed=transposed)
})
