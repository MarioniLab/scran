#' Approximate silhouette width
#'
#' \linkS4class{SingleCellExperiment}-compatible wrapper around \code{\link{approxSilhouette}},
#' for quickly computing an approximate silhouette width for each cell.
#'
#' @inheritParams clusterPurity
#' 
#' @return
#' A \linkS4class{DataFrame} of purity statistics where each row corresponds to a cell in \code{x}, 
#' see \code{?\link{approxSilhouette}} for details.
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
#' \code{\link{approxSilhouette}}, for the details of the approximate calculations.
#'
#' \code{\link{clusterPurity}}, for another method of evaluating cluster separation.
#'
#' @name clusterSilhouette
NULL

#' @importFrom Matrix t
#' @importFrom bluster approxSilhouette
.cluster_silhouette <- function(x, ..., transposed=FALSE, subset.row=NULL) {
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    }
    approxSilhouette(x, ...)
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
