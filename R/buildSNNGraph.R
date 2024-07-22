#' Build a nearest-neighbor graph
#'
#' \linkS4class{SingleCellExperiment}-friendly wrapper around the \code{\link{makeSNNGraph}} and \code{\link{makeKNNGraph}} functions for creating nearest-neighbor graphs.
#'
#' @param x A matrix-like object containing expression values for each gene (row) in each cell (column).
#' These dimensions can be transposed if \code{transposed=TRUE}.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such an expression matrix.
#' If \code{x} is a SingleCellExperiment and \code{use.dimred} is set, its \code{\link{reducedDims}} will be used instead.
#' @param d An integer scalar specifying the number of dimensions to use for a PCA on the expression matrix prior to the nearest neighbor search.
#' Ignored for the ANY method if \code{transposed=TRUE} and for the SingleCellExperiment methods if \code{use.dimred} is set.
#' @param transposed A logical scalar indicating whether \code{x} is transposed (i.e., rows are cells).
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' Only used when \code{transposed=FALSE}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, if \code{d} is not \code{NA}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object to use for parallel processing.
#' @param ... For the generics, additional arguments to pass to the specific methods.
#' 
#' For the ANY methods, additional arguments to pass to \code{\link{makeSNNGraph}} or \code{\link{makeKNNGraph}}.
#'
#' For the SummarizedExperiment methods, additional arguments to pass to the corresponding ANY method.
#'
#' For the SingleCellExperiment methods, additional arguments to pass to the corresponding SummarizedExperiment method.
#' @param assay.type A string specifying which assay values to use.
#' @param use.dimred A string specifying whether existing values in \code{reducedDims(x)} should be used.
#' 
#' @return
#' A \link{graph} where nodes are cells and edges represent connections between nearest neighbors,
#' see \code{?\link{makeSNNGraph}} for more details.
#' 
#' @author 
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{makeSNNGraph}} and \code{\link{makeKNNGraph}}, for the underlying functions that do the work.
#' 
#' See \code{\link{cluster_walktrap}} and related functions in \pkg{igraph} for clustering based on the produced graph.
#'
#' \code{\link{clusterCells}}, for a more succinct way of performing graph-based clustering.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ncells=500)
#' sce <- logNormCounts(sce)
#'
#' g <- buildSNNGraph(sce)
#' clusters <- igraph::cluster_fast_greedy(g)$membership
#' table(clusters)
#'
#' # Any clustering method from igraph can be used:
#' clusters <- igraph::cluster_walktrap(g)$membership
#' table(clusters)
#'
#' # Smaller 'k' usually yields finer clusters:
#' g <- buildSNNGraph(sce, k=5)
#' clusters <- igraph::cluster_walktrap(g)$membership
#' table(clusters)
#'
#' # Graph can be built off existing reducedDims results:
#' sce <- scater::runPCA(sce)
#' g <- buildSNNGraph(sce, use.dimred="PCA")
#' clusters <- igraph::cluster_fast_greedy(g)$membership
#' table(clusters)
#' 
#' @name buildSNNGraph
#' @docType methods
NULL

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular bsparam
#' @importFrom bluster makeSNNGraph
.buildSNNGraph <- function(x, ..., d=50, transposed=FALSE, subset.row=NULL, BSPARAM=bsparam(), BPPARAM=SerialParam()) { 
    input <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
    makeSNNGraph(input, ..., BPPARAM=BPPARAM)
}

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular bsparam
#' @importFrom bluster makeKNNGraph
.buildKNNGraph <- function(x, ..., d=50, transposed=FALSE, subset.row=NULL, BSPARAM=bsparam(), BPPARAM=SerialParam()) { 
    input <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
    makeKNNGraph(input, ..., BPPARAM=BPPARAM)
}

.setup_knn_data <- function(x, subset.row, d, transposed, k, BNPARAM, BSPARAM, BPPARAM) {
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    } 
    
    # Reducing dimensions, if 'd' is less than the number of genes.
    if (!is.na(d) && d < ncol(x)) {
        svd.out <- .centered_SVD(x, max.rank=d, keep.right=FALSE, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        x <- .svd_to_pca(svd.out, d, named=FALSE)
    }

    x
}

#########################
# S4 method definitions #
#########################

#' @export
#' @rdname buildSNNGraph
setGeneric("buildSNNGraph", function(x, ...) standardGeneric("buildSNNGraph"))

#' @export
#' @rdname buildSNNGraph
setMethod("buildSNNGraph", "ANY", .buildSNNGraph)

#' @export
#' @rdname buildSNNGraph
#' @importFrom SummarizedExperiment assay
setMethod("buildSNNGraph", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .buildSNNGraph(assay(x, i=assay.type), transposed=FALSE, ...)
})

#' @export
#' @rdname buildSNNGraph
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim 
setMethod("buildSNNGraph", "SingleCellExperiment", function(x, ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .buildSNNGraph(reducedDim(x, use.dimred), d=NA, transposed=TRUE, ...)
    } else {
        callNextMethod(x=x, ...)
    }
})

#' @export
#' @rdname buildSNNGraph
setGeneric("buildKNNGraph", function(x, ...) standardGeneric("buildKNNGraph"))

#' @export
#' @rdname buildSNNGraph
setMethod("buildKNNGraph", "ANY", .buildKNNGraph)

#' @export
#' @rdname buildSNNGraph
#' @importFrom SummarizedExperiment assay
setMethod("buildKNNGraph", "SingleCellExperiment", function(x, ..., assay.type="logcounts"){ 
    .buildKNNGraph(assay(x, i=assay.type), transposed=FALSE, ...)
})

#' @export
#' @rdname buildSNNGraph
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim 
setMethod("buildKNNGraph", "SingleCellExperiment", function(x, ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .buildKNNGraph(reducedDim(x, use.dimred), d=NA, transposed=TRUE, ...)
    } else {
        callNextMethod(x=x, ...)
    }
})
