#' Wrappers for graph-based clustering
#'
#' Perform graph-based clustering using community detection methods on a nearest-neighbor graph,
#' where nodes represent cells or k-means centroids.
#'
#' @inheritParams buildSNNGraph
#' @param use.kmeans Logical scalar indicating whether k-means clustering should be performed.
#' @param clusterFUN Function that can take a \link{graph} object and return a \link{communities} object,
#' see examples in the \pkg{igraph} package.
#' @param kmeans.centers Integer scalar specifying the number of clusters to use for k-means clustering.
#' Defaults to the square root of the number of cells in \code{x}.
#' @param kmeans.args List containing additional named arguments to pass to \code{\link{kmeans}}.
#' @param full.stats Logical scalar indicating whether to return more statistics regarding the k-means clustering. 
#' @param ... For the generics, additional arguments to pass to the specific methods.
#'
#' For the ANY methods, additional arguments to pass to \code{\link{buildSNNGraph}} or \code{\link{buildKNNGraph}}.
#'
#' For the SummarizedExperiment method, additional arguments to pass to the corresponding ANY method.
#'
#' For the SingleCellExperiment method, additional arguments to pass to the corresponding SummarizedExperiment method.
#'
#' @return 
#' If \code{full.stats=FALSE}, a factor is returned containing the cluster assignment for each cell.
#' 
#' If \code{full.stats=TRUE} and \code{use.kmeans=TRUE}, a \linkS4class{DataFrame} is returned with one row per cell.
#' This contains the columns \code{kmeans}, specifying the assignment of each cell to a k-means cluster;
#' and \code{igraph}, specifying the assignment of each cell to a graph-based cluster operating on the k-means clusters.
#' In addition, the \code{\link{metadata}} contains \code{graph}, a \link{graph} object where each node is a k-means cluster;
#' and \code{membership}, the graph-based cluster to which each node is assigned.
#'
#' @details
#' By default, these functions simply call \code{\link{buildSNNGraph}} or \code{\link{buildKNNGraph}}
#' followed by the specified \code{clusterFUN} to generate a clustering.
#' We use the Walktrap algorithm by default as it has a cool-sounding name,
#' but users can feel free to swap it for some other algorithm (e.g., \code{\link{cluster_louvain}}).
#'
#' For large datasets, we can perform vector quantization with k-means
#' to create centroids that are then subjected to graph-based clustering.
#' The label for each cell is then defined as the label of the centroid to which it was assigned.
#' In this approach, k-means and graph-based clustering complement each other;
#' the former improves computational efficiency and mitigates density-dependent dilation,
#' while the latter aggregates the centroids for easier interpretation.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{buildSNNGraph}} and \code{\link{buildKNNGraph}}, to build the nearest-neighbor graphs.
#'
#' \code{\link{cluster_walktrap}} and related functions, to detect communities within those graphs.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ncells=500)
#' sce <- logNormCounts(sce)
#'
#' clusters <- clusterSNNGraph(sce)
#' table(clusters)
#'
#' # Can pass usual arguments to buildSNNGraph:
#' clusters2 <- clusterSNNGraph(sce, k=5)
#' table(clusters2)
#'
#' # Works with low-dimensional inputs:
#' sce <- scater::runPCA(sce, ncomponents=10)
#' clusters3 <- clusterSNNGraph(sce, use.dimred="PCA")
#' table(clusters3)
#'
#' # Turn on k-means for larger datasets, e.g., 
#' # assuming we already have a PCA result:
#' set.seed(101010)
#' bigpc <- matrix(rnorm(2000000), ncol=20)
#' clusters4 <- clusterSNNGraph(bigpc, d=NA, use.kmeans=TRUE, transposed=TRUE)
#' table(clusters4)
#'
#' # Extract the graph for more details:
#' clusters5 <- clusterSNNGraph(sce, use.dimred="PCA", 
#'     use.kmeans=TRUE, full.stats=TRUE)
#' head(clusters5)
#' metadata(clusters5)$graph
#'
#' @name clusterSNNGraph
NULL

#' @importFrom Matrix t
#' @importFrom stats kmeans 
#' @importFrom igraph V V<-
.clusterXNNGraph <- function(x, ..., subset.row=NULL, transposed=FALSE, use.kmeans=FALSE, kmeans.centers=NULL, 
    kmeans.args=list(), full.stats=FALSE, graphFUN, clusterFUN) 
{
    if (use.kmeans) {
        if (is.null(kmeans.centers)) {
            kmeans.centers <- ceiling(sqrt(nrow(x)))
        }

        if (!transposed) {
            if (!is.null(subset.row)) {
                x <- x[subset.row,,drop=FALSE]
            }
            x <- t(x)
        }

        kmeans.args <- c(list(as.matrix(x), centers=kmeans.centers), kmeans.args)
        k.out <- do.call(kmeans, kmeans.args)

        nn.g <- graphFUN(k.out$centers, ..., transposed=TRUE)
        clust.g <- clusterFUN(nn.g)$membership
        output <- clust.g[k.out$cluster]

        if (full.stats) {
            df <- DataFrame(
                kmeans=k.out$cluster,
                igraph=factor(output)
            )
            metadata(df) <- list(graph=nn.g, membership=clust.g)
            return(df)
        }

    } else {
        nn.g <- graphFUN(x, ..., subset.row=subset.row, transposed=transposed)
        output <- clusterFUN(nn.g)$membership
    }

    factor(output)
}

#' @importFrom igraph cluster_walktrap
.clusterKNNGraph <- function(x, ..., clusterFUN=cluster_walktrap, subset.row=NULL, transposed=FALSE, 
    use.kmeans=FALSE, kmeans.centers=NULL, kmeans.args=list(), full.stats=FALSE) 
{
    .clusterXNNGraph(x=x, ..., clusterFUN=clusterFUN, subset.row=subset.row, transposed=transposed,
        use.kmeans=use.kmeans, kmeans.centers=kmeans.centers, kmeans.args=kmeans.args,
        full.stats=full.stats, graphFUN=buildKNNGraph)
}

#' @importFrom igraph cluster_walktrap
.clusterSNNGraph <- function(x, ..., clusterFUN=cluster_walktrap, subset.row=NULL, transposed=FALSE, 
    use.kmeans=FALSE, kmeans.centers=NULL, kmeans.args=list(), full.stats=FALSE) 
{
    .clusterXNNGraph(x=x, ..., clusterFUN=clusterFUN, subset.row=subset.row, transposed=transposed,
        use.kmeans=use.kmeans, kmeans.centers=kmeans.centers, kmeans.args=kmeans.args,
        full.stats=full.stats, graphFUN=buildSNNGraph)
}


#########################
# S4 method definitions #
#########################

#' @export
#' @rdname clusterSNNGraph
setGeneric("clusterSNNGraph", function(x, ...) standardGeneric("clusterSNNGraph"))

#' @export
#' @rdname clusterSNNGraph
setMethod("clusterSNNGraph", "ANY", .clusterSNNGraph)

#' @export
#' @rdname clusterSNNGraph
#' @importFrom SummarizedExperiment assay
setMethod("clusterSNNGraph", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .clusterSNNGraph(assay(x, assay.type), ..., transposed=FALSE)
})

#' @export
#' @rdname clusterSNNGraph
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
setMethod("clusterSNNGraph", "SingleCellExperiment", function(x, ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .clusterSNNGraph(reducedDim(x, use.dimred), ..., d=NA, transposed=TRUE)
    } else {
        callNextMethod(x=x, ...)
    }
})

#' @export
#' @rdname clusterSNNGraph
setGeneric("clusterKNNGraph", function(x, ...) standardGeneric("clusterKNNGraph"))

#' @export
#' @rdname clusterSNNGraph
setMethod("clusterKNNGraph", "ANY", .clusterKNNGraph)

#' @export
#' @rdname clusterSNNGraph
#' @importFrom SummarizedExperiment assay
setMethod("clusterKNNGraph", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .clusterKNNGraph(assay(x, assay.type), ..., transposed=FALSE)
})

#' @export
#' @rdname clusterSNNGraph
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
setMethod("clusterKNNGraph", "SingleCellExperiment", function(x, ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .clusterKNNGraph(reducedDim(x, use.dimred), ..., d=NA, transposed=TRUE)
    } else {
        callNextMethod(x=x, ...)
    }
})
