#' Quick MST-based pseudotime
#'
#' A convenience wrapper to quickly compute a minimum spanning tree (MST) on the cluster centroids
#' to obtain a pseudotime ordering of the cells.
#'
#' @param x A named list of numeric matrices containing dimensionality reduction results.
#' All matrices should have the same number of cells, i.e., rows.
#' Alternatively, a \linkS4class{SingleCellExperiment} containing such results in its \code{\link{reducedDims}}.
#' @param clusters A factor of length equal to the number of cells in \code{x},
#' specifying the cluster assignment for each cell.
#' @param use Integer scalar or string specifying the entry of \code{x} to use for MST construction and pseudotime calculations.
#' @param start Arguments passed to \code{\link{orderClusterMST}}.
#' @param outgroup,outscale Arguments passed to \code{\link{createClusterMST}}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SingleCellExperiment method, further arguments to pass to the ANY method.
#'
#' @details
#' This function simply calls, in order:
#' \itemize{
#' \item \code{\link{rowsum}}, to compute the average low-dimensional coordinates for each cluster.
#' \item \code{\link{createClusterMST}} on the average coordinates specified by \code{use}.
#' \item \code{\link{orderClusterMST}} on the average and per-cell coordinates specified by \code{use}.
#' \item \code{\link{connectClusterMST}} on the average coordinates for all entries of \code{x}.
#' }
#'
#' @return
#' A \linkS4class{List} containing:
#' \itemize{
#' \item \code{centered}, a list of numeric matrices containing the averaged coordinates for each cluster.
#' Each matrix corresponds to a dimensionality reduction result in \code{x}.
#' \item \code{mst}, a \link{graph} object containing the cluster-level MST computed on the coordinates from \code{use}.
#' \item \code{ordering}, a numeric matrix of pseudotimes for various paths through the MST computed from \code{use}.
#' \item \code{connected}, a list of data.frames containing the edge coordinates between centers.
#' Each data.frame corresponds to a dimensionality reduction result in \code{x}.
#' }
#'
#' @seealso
#' \code{\link{createClusterMST}} and friends, for the functions that do the actual work.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up some data:
#' library(scuttle)
#' sce <- mockSCE(ncells=500)
#' sce <- logNormCounts(sce)
#' sce <- scater::runPCA(sce)
#' clusters <- clusterSNNGraph(sce, use.dimred="PCA")
#'
#' # Quickly computing the pseudotime:
#' out <- quickPseudotime(sce, clusters, use="PCA")
#' out$mst
#' head(out$ordering)
#'
#' @name quickPseudotime
#' @aliases
#' quickPseudotime,ANY-method
#' quickPseudotime,SingleCellExperiment-method 
NULL

#' @importFrom SingleCellExperiment reducedDims reducedDims<-
#' @importFrom SingleCellExperiment SingleCellExperiment
.quick_pseudotime <- function(x, clusters, use=1, outgroup=FALSE, outscale=3, start=NULL) {
    tab <- table(clusters)
    centered <- x
    for (i in seq_along(x)) {
        current <- rowsum(x[[i]], clusters)
        centered[[i]] <- current/as.integer(tab[rownames(current)])
    }

    mst <- createClusterMST(centered[[use]], outgroup=outgroup, outscale=outscale)
    connected <- lapply(centered, FUN=connectClusterMST, mst=mst)
    ordering <- orderClusterMST(x[[use]], ids=clusters, centers=centered[[use]], mst=mst, start=start)

    List(
        centered=centered,
        mst=mst,
        ordering=ordering,
        connected=connected
    )
}

#' @export
#' @rdname quickPseudotime
setGeneric("quickPseudotime", function(x, ...) standardGeneric("quickPseudotime"))

#' @export
#' @rdname quickPseudotime
setMethod("quickPseudotime", "ANY", .quick_pseudotime)

#' @export
#' @rdname quickPseudotime
#' @importFrom SingleCellExperiment colLabels reducedDims
setMethod("quickPseudotime", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), ...) {
    .quick_pseudotime(reducedDims(x), clusters=clusters, ...)
})
