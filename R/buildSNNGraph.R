#' Build a nearest-neighbor graph
#'
#' Build a shared or k-nearest-neighbors graph of cells based on similarities in their expression profiles.
#'
#' @param x A matrix-like object containing expression values for each gene (row) in each cell (column).
#' These dimensions can be transposed if \code{transposed=TRUE}.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such an expression matrix.
#'
#' If \code{x} is a SingleCellExperiment and \code{use.dimred} is set,
#' graph building will be performed from its \code{\link{reducedDims}}.
#' @param k An integer scalar specifying the number of nearest neighbors to consider during graph construction.
#' @param d An integer scalar specifying the number of dimensions to use for the search.
#' Ignored for the SingleCellExperiment methods if \code{use.dimred} is set.
#' @param type A string specifying the type of weighting scheme to use for shared neighbors.
#' @param directed A logical scalar indicating whether the output of \code{buildKNNGraph} should be a directed graph.
#' @param transposed A logical scalar indicating whether \code{x} is transposed (i.e., rows are cells).
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' Only used when \code{transposed=FALSE}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, if \code{d} is not \code{NA}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object to use for parallel processing.
#' @param ... For the generics, additional arguments to pass to the specific methods.
#' 
#' For the SummarizedExperiment methods, additional arguments to pass to the corresponding ANY method.
#'
#' For the SingleCellExperiment methods, additional arguments to pass to the corresponding SummarizedExperiment method.
#' @param assay.type A string specifying which assay values to use.
#' @param use.dimred A string specifying whether existing values in \code{reducedDims(x)} should be used.
#' @param indices An integer matrix where each row corresponds to a cell and contains the indices of the \code{k} nearest neighbors (by increasing distance) from that cell.
#' 
#' @details
#' The \code{buildSNNGraph} method builds a shared nearest-neighbour graph using cells as nodes.
#' For each cell, its \code{k} nearest neighbours are identified using the \code{\link{findKNN}} function,
#' based on distances between their expression profiles (Euclidean by default).
#' An edge is drawn between all pairs of cells that share at least one neighbour,
#' weighted by the characteristics of the shared nearest neighbors:
#' \itemize{
#' \item If \code{type="rank"}, the weighting scheme defined by Xu and Su (2015) is used.
#' The weight between two nodes is \eqn{k - r/2} where \eqn{r} is the smallest sum of ranks for any shared neighboring node.
#' For example, if one node was the closest neighbor of each of two nodes, the weight between the two latter nodes would be \eqn{k - 1}.
#' For the purposes of this ranking, each node has a rank of zero in its own nearest-neighbor set. 
#' More shared neighbors, or shared neighbors that are close to both cells, will generally yield larger weights.
#' \item If \code{type="number"}, the weight between two nodes is simply the number of shared nearest neighbors between them.
#' The weight can range from zero to \eqn{k + 1}, as the node itself is included in its own nearest-neighbor set.
#' This is a simpler scheme that is also slightly faster but does not account for the ranking of neighbors within each set.
#' \item If \code{type="jaccard"}, the weight between two nodes is the Jaccard similarity between the two sets of neighbors for those nodes.
#' This weight can range from zero to 1, and is a monotonic transformation of the weight used by \code{type="number"}.
#' It is provided for consistency with other clustering algorithms such as those in \pkg{seurat}.
#' }
#' 
#' The aim is to use the SNN graph to perform clustering of cells via community detection algorithms in the \pkg{igraph} package.
#' This is faster and more memory efficient than hierarchical clustering for large numbers of cells.
#' In particular, it avoids the need to construct a distance matrix for all pairs of cells.
#' Only the identities of nearest neighbours are required, which can be obtained quickly with methods in the \pkg{BiocNeighbors} package.
#' 
#' The choice of \code{k} controls the connectivity of the graph and the resolution of community detection algorithms.
#' Smaller values of \code{k} will generally yield smaller, finer clusters, while increasing \code{k} will increase the connectivity of the graph and make it more difficult to resolve different communities.
#' The value of \code{k} can be roughly interpreted as the anticipated size of the smallest subpopulation.
#' If a subpopulation in the data has fewer than \code{k+1} cells, \code{buildSNNGraph} and \code{buildKNNGraph} will forcibly construct edges between cells in that subpopulation and cells in other subpopulations. 
#' This increases the risk that the subpopulation will not form its own cluster as it is more interconnected with the rest of the cells in the dataset.
#' 
#' Note that the setting of \code{k} here is slightly different from that used in SNN-Cliq.
#' The original implementation considers each cell to be its first nearest neighbor that contributes to \code{k}.
#' In \code{buildSNNGraph}, the \code{k} nearest neighbours refers to the number of \emph{other} cells.
#' 
#' The \code{buildKNNGraph} method builds a simpler k-nearest neighbour graph.
#' Cells are again nodes, and edges are drawn between each cell and its k-nearest neighbours.
#' No weighting of the edges is performed.
#' In theory, these graphs are directed as nearest neighour relationships may not be reciprocal.
#' However, by default, \code{directed=FALSE} such that an undirected graph is returned.
#' 
#' @section Choice of input data:
#' In practice, PCA is performed on a matrix-like \code{x} to obtain the first \code{d} principal components.
#' The k-NN search can then be rapidly performed on the PCs rather than on the full expression matrix.
#' By default, the first 50 components are chosen, which should retain most of the substructure in the data set.
#' If \code{d} is \code{NA} or greater than or equal to the number of cells, no dimensionality reduction is performed.
#'     
#' The PCA is performed using methods the \code{\link{runSVD}} function from the \pkg{BiocSingular} package.
#' To improve speed, this can be done using approximate algorithms by modifying \code{BSPARAM}, e.g., to \code{\link{IrlbaParam}()}.
#' Approximate algorithms will converge towards the correct result but often involve some random initialization and thus are technically dependent on the session seed.
#' For full reproducibility, users are advised to call \code{\link{set.seed}} beforehand when using this option.
#' 
#' Expression values in \code{x} should typically be on the log-scale, e.g., log-transformed counts.
#' Ranks can also be used for greater robustness, e.g., from \code{\link{scaledColRanks}}. 
#' (Dimensionality reduction is still okay when ranks are provided - running PCA on ranks is equivalent to running MDS on the distance matrix derived from Spearman's rho.)
#' 
#' If the input matrix \code{x} is already transposed for the ANY method, \code{transposed=TRUE} avoids an unnecessary internal transposition.
#' A typical use case is when \code{x} contains some reduced dimension coordinates with cells in the rows.
#' In such cases, setting \code{transposed=TRUE} and \code{d=NA} will use the input coordinates directly for graph-building.
#' 
#' The same principles apply when \code{x} is a \linkS4class{SingleCellExperiment} object, 
#' except that the relevant matrix is now retrieved from the assays using \code{assay.type}. 
#' If \code{use.dimred} is not \code{NULL}, existing PCs are used from the specified entry of \code{reducedDims(x)}, 
#' and any setting of \code{d} and \code{subset.row} are ignored.
#'
#' The \code{neighborsToSNNGraph} and \code{neighborsToKNNGraph} functions operate directly on a matrix of nearest neighbor indices,
#' obtained using functions like \code{\link{findKNN}}.
#' This may be useful for constructing a graph from precomputed nearest-neighbor search results.
#' Note that the user is responsible for ensuring that the indices are valid (i.e., \code{range(indices)} is positive and no greater than \code{max(indices)}).
#' 
#' @return
#' A \link{graph} where nodes are cells and edges represent connections between nearest neighbors.
#' For \code{buildSNNGraph}, these edges are weighted by the number of shared nearest neighbors.
#' For \code{buildKNNGraph}, edges are not weighted but may be directed if \code{directed=TRUE}.
#' 
#' @author 
#' Aaron Lun, with KNN code contributed by Jonathan Griffiths.
#' 
#' @seealso
#' See \code{\link{make_graph}} for details on the graph output object.
#' 
#' See \code{\link{cluster_walktrap}}, \code{\link{cluster_louvain}} and related functions in \pkg{igraph} for clustering based on the produced graph.
#' 
#' Also see \code{\link{findKNN}} for specifics of the nearest-neighbor search.
#' 
#' @references
#' Xu C and Su Z (2015).
#' Identification of cell types from single-cell transcriptomes using a novel clustering method.
#' \emph{Bioinformatics} 31:1974-80
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
NULL

#' @importFrom BiocNeighbors KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular bsparam
.buildSNNGraph <- function(x, k=10, d=50, 
    type=c("rank", "number", "jaccard"),
    transposed=FALSE, subset.row=NULL, 
    BNPARAM=KmknnParam(), BSPARAM=bsparam(), BPPARAM=SerialParam()) 
{ 
    nn.out <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed,
        k=k, BNPARAM=BNPARAM, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
    neighborsToSNNGraph(nn.out$index, type=match.arg(type))
}

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular bsparam
#' @importFrom BiocNeighbors KmknnParam
.buildKNNGraph <- function(x, k=10, d=50, directed=FALSE, transposed=FALSE,
    subset.row=NULL, BNPARAM=KmknnParam(), BSPARAM=bsparam(), BPPARAM=SerialParam()) 
{ 
    nn.out <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed,
        k=k, BNPARAM=BNPARAM, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
    neighborsToKNNGraph(nn.out$index, directed=directed)
}

#' @importFrom BiocNeighbors findKNN
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
   
    # Finding the KNNs. 
    findKNN(x, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=FALSE)
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

############################
# Graph-building functions #
############################

#' @export
#' @rdname buildSNNGraph
#' @importFrom igraph make_graph simplify "E<-"
neighborsToSNNGraph <- function(indices, type=c("rank", "number", "jaccard")) {
    type <- match.arg(type)
    if (type=="rank") {
        g.out <- build_snn_rank(indices)
    } else {
        g.out <- build_snn_number(indices)
    }

    edges <- g.out[[1]] 
    weights <- g.out[[2]]
    if (type=="jaccard") {
        weights <- weights / (2 * (ncol(indices) + 1) - weights)
    }

    g <- make_graph(edges, directed=FALSE)
    E(g)$weight <- weights
    simplify(g, edge.attr.comb="first") # symmetric, so doesn't really matter.
}


#' @export
#' @rdname buildSNNGraph
#' @importFrom igraph make_graph simplify
neighborsToKNNGraph <- function(indices, directed=FALSE) {
    start <- as.vector(row(indices))
    end <- as.vector(indices)
    interleaved <- as.vector(rbind(start, end))
    
    if (directed) { 
        g <- make_graph(interleaved, directed=TRUE)
    } else {
        g <- make_graph(interleaved, directed=FALSE)
        g <- simplify(g, edge.attr.comb = "first")
    }
    g
}
