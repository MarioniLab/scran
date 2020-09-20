#' Evaluate cluster purity
#'
#' Determine whether cells are surrounded by neighbors that are assigned to the same cluster.
#' This function has now been deprecated in favor of \code{\link{neighborPurity}} from the \pkg{bluster} package.
#'
#' @inheritParams buildSNNGraph
#' @param clusters A vector or factor of cluster IDs to pass to \code{\link{neighborPurity}}.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the ANY method, arguments to pass to \code{\link{neighborPurity}}.
#'
#' For the SummarizedExperiment method, arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, arguments to pass to the SummarizedExperiment method.
#'
#' @return
#' A \linkS4class{DataFrame} of purity statistics where each row corresponds to a cell in \code{x}, 
#' see \code{?\link{neighborPurity}} for details.
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
#' @seealso
#' \code{\link{approxSilhouette}}, for another method of evaluating cluster separation.
#'
#' @name clusterPurity
NULL

#' @importFrom Matrix t
#' @importFrom bluster neighborPurity
.cluster_purity <- function(x, ..., transposed=FALSE, subset.row=NULL) {
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    }
    .Deprecated(old="clusterPurity", new="bluster::neighborPurity")
    neighborPurity(x, ...)
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
