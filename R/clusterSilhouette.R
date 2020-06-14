#' Approximate silhouette width
#'
#' Given a clustering, compute a fast approximate silhouette width for each cell.
#'
#' @inheritParams clusterPurity
#' 
#' @details
#' The silhouette width is a general purpose method for evaluating the separation between clusters 
#' but requires calculating the average distance between pairs of cells within or between clusters.
#' This function approximates the average distance with \eqn{d}, defined as the larger of:
#' \itemize{
#' \item The distance from the current cell to the centroid of cluster \eqn{X}.
#' This is a reasonable approximation when the cell is distant to \eqn{X} relative to the latter's variation.
#' \item The root mean squared distance between cells in cluster \eqn{X}.
#' This is a reasonable approximation when the cell lies inside the bulk of points for \eqn{X}.
#' }
#' The approximate silhouette width for each cell can then be calculated with two values of \eqn{d},
#' computed by setting \eqn{X} to the cluster of the current cell or the closest other cluster.
#' 
#' @return
#' A \linkS4class{DataFrame} with one row per cell in \code{x} and the columns:
#' \itemize{
#' \item \code{width}, a numeric field containing the approximate silhouette width of the current cell.
#' \item \code{other}, the closest cluster other than the one to which the current cell is assigned.
#' }
#' Row names are defined as the column names of \code{x}.
#' 
#' @author Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' g <- buildSNNGraph(sce)
#' clusters <- igraph::cluster_walktrap(g)$membership
#' out <- clusterSilhouette(sce, clusters)
#' boxplot(split(out$width, clusters))
#'
#' # Mocking up a stronger example:
#' ngenes <- 1000
#' centers <- matrix(rnorm(ngenes*3), ncol=3)
#' clusters <- sample(1:3, ncol(sce), replace=TRUE)
#'
#' y <- centers[,clusters]
#' y <- y + rnorm(length(y))
#' 
#' out2 <- clusterSilhouette(y, clusters)
#' boxplot(split(out2$width, clusters))
#'
#' @seealso
#' \code{silhouette} from the \pkg{cluster} package, for the exact calculation.
#'
#' \code{\link{clusterPurity}}, for another method of evaluating cluster separation.
#'
#' @name clusterSilhouette
NULL

#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix t
#' @importFrom BiocNeighbors queryKNN
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedMatrixStats colVars
.cluster_silhouette <- function(x, clusters, k=50, transposed=FALSE, subset.row=NULL, BPPARAM=SerialParam()) {
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    }

    x <- as.matrix(x)
    uclust <- sort(unique(clusters))
    averaged <- matrix(0, length(uclust), ncol(x))
    self.dist <- numeric(nrow(x))
    clust.rmsd <- numeric(length(uclust))

    for (i in seq_along(uclust)) {
        current <- uclust[i]==clusters
        xcurrent <- x[current,,drop=FALSE]
        centroid <- colMeans(xcurrent)
        averaged[i,] <- centroid

        # Technically, this is the root of the estimate of the expectation of
        # the squared distances. Recall that we have (X - Y)^2 for i.i.d. RV's 
        # X and Y; this expands to (X - u)^2 - 2*(X - u)(Y - u) + (Y - u)^2,
        # the expectation of which is simply twice the variance of X or Y. 
        cur.rmsd <- sqrt(sum(colVars(xcurrent, center=centroid)) * 2)
        clust.rmsd[i] <- cur.rmsd

        delta2center <- sweep(xcurrent, 2, centroid, FUN='-')
        dist2center <- sqrt(rowSums(delta2center^2))
        self.dist[current] <- pmax(dist2center, cur.rmsd)
    }

    # Getting the next-closest (grabbing 2 and ruling out the self).
    closest <- queryKNN(query=x, X=averaged, k=2, BPPARAM=BPPARAM)
    not.self <- uclust[closest$index[,1]]!=clusters
    other.clust <- ifelse(not.self, closest$index[,1], closest$index[,2])
    other.dist <- ifelse(not.self, closest$distance[,1], closest$distance[,2])
    other.dist <- pmax(other.dist, clust.rmsd[other.clust])

    DataFrame(
        width=(other.dist - self.dist)/pmax(other.dist, self.dist),
        other=uclust[other.clust],
        row.names=rownames(x)
    )
}

#' @export
#' @rdname clusterSilhouette
setGeneric("clusterSilhouette", function(x, ...) standardGeneric("clusterSilhouette"))

#' @export
#' @rdname clusterSilhouette
setMethod("clusterSilhouette", "ANY", .cluster_silhouette)

#' @export
#' @rdname clusterSilhouette
#' @importFrom SummarizedExperiment assay
setMethod("clusterSilhouette", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .cluster_silhouette(assay(x, assay.type), ..., transposed=FALSE)
})

#' @export
#' @rdname clusterSilhouette
#' @importFrom SingleCellExperiment reducedDim colLabels
#' @importFrom SummarizedExperiment assay
setMethod("clusterSilhouette", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"),
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
    .cluster_silhouette(x, clusters=clusters, ..., transposed=transposed)
})
