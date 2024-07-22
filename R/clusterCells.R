#' Cluster cells in a SingleCellExperiment
#'
#' A \linkS4class{SingleCellExperiment}-compatible wrapper around \code{\link{clusterRows}} from the \pkg{bluster} package.
#'
#' @param x A \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing cells in the columns.
#' @param assay.type Integer or string specifying the assay values to use for clustering, typically log-normalized expression.
#' @param use.dimred Integer or string specifying the reduced dimensions to use for clustering, typically PC scores.
#' Only used when \code{assay.type=NULL}, and only applicable if \code{x} is a SingleCellExperiment.
#' @param BLUSPARAM A \linkS4class{BlusterParam} object specifying the clustering algorithm to use,
#' defaults to a graph-based method.
#' @param ... Further arguments to pass to \code{\link{clusterRows}}.
#'
#' @return A factor of cluster identities for each cell in \code{x},
#' or a list containing such a factor - see the return value of \code{?\link{clusterRows}}.
#'
#' @details
#' This is largely a convenience wrapper to avoid the need to manually extract the relevant assays or reduced dimensions from \code{x}.
#' Altering \code{BLUSPARAM} can easily change the parameters or algorithm used for clustering -
#' see \code{?"\link{BlusterParam-class}"} for more details.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' # From log-expression values:
#' clusters <- clusterCells(sce, assay.type="logcounts")
#'
#' # From PCs:
#' sce <- scater::runPCA(sce)
#' clusters2 <- clusterCells(sce, use.dimred="PCA")
#'
#' # With different parameters:
#' library(bluster)
#' clusters3 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=NNGraphParam(k=5))
#' 
#' # With different algorithms:
#' clusters4 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=KmeansParam(centers=10))
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
#' @importFrom bluster clusterRows NNGraphParam
clusterCells <- function(x, assay.type=NULL, use.dimred=NULL, BLUSPARAM=NNGraphParam(), ...) {
    if (!is.null(assay.type)) {
        x <- t(assay(x, assay.type))
    } else if (!is.null(use.dimred)) {
        x <- reducedDim(x, use.dimred)
    } else {
        stop("either 'assay.type=' or 'use.dimred=' must be specified")
    }

    x <- as.matrix(x)
    clusterRows(x, BLUSPARAM=BLUSPARAM, ...)
}
