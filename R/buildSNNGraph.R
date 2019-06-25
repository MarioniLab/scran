#' @importFrom igraph make_graph simplify "E<-"
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular ExactParam
.buildSNNGraph <- function(x, k=10, d=50, 
    type=c("rank", "number", "jaccard"),
    transposed=FALSE, subset.row=NULL, 
    BNPARAM=KmknnParam(), BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Builds a shared nearest-neighbor graph, where edges are present between each 
# cell and any other cell with which it shares at least one neighbour. Each edge
# is weighted based on the ranks of the shared nearest neighbours of the two cells, 
# as described in the SNN-Cliq paper; or by the number of shared numbers;
# or by the Jaccard index, a la Seurat's default.
#
# written by Aaron Lun
# created 3 April 2017
{ 
    nn.out <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed,
        k=k, BNPARAM=BNPARAM, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 

    # Building the SNN graph.
    type <- match.arg(type)
    if (type=="rank") {
        g.out <- .Call(cxx_build_snn_rank, nn.out$index)
    } else {
        g.out <- .Call(cxx_build_snn_number, nn.out$index)
    }

    edges <- g.out[[1]] 
    weights <- g.out[[2]]
    if (type=="jaccard") {
        weights <- weights / (2 * (k + 1) - weights)
    }

    g <- make_graph(edges, directed=FALSE)
    E(g)$weight <- weights
    g <- simplify(g, edge.attr.comb="first") # symmetric, so doesn't really matter.
    return(g)
}

#' @importFrom igraph make_graph simplify
#' @importFrom BiocParallel SerialParam
.buildKNNGraph <- function(x, k=10, d=50, directed=FALSE, transposed=FALSE,
    subset.row=NULL, BNPARAM=KmknnParam(), BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Builds a k-nearest-neighbour graph, where edges are present between each
# cell and its 'k' nearest neighbours. Undirected unless specified otherwise.
#
# written by Aaron Lun, Jonathan Griffiths
# created 16 November 2017
{ 
    nn.out <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed,
        k=k, BNPARAM=BNPARAM, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 

    # Building the KNN graph.
    start <- as.vector(row(nn.out$index))
    end <- as.vector(nn.out$index)
    interleaved <- as.vector(rbind(start, end))
    
    if (directed) { 
        g <- make_graph(interleaved, directed=TRUE)
    } else {
        g <- make_graph(interleaved, directed=FALSE)
        g <- simplify(g, edge.attr.comb = "first")
    }
    return(g)
}

######################
# Internal functions #
######################

#' @importFrom stats prcomp 
#' @importFrom BiocNeighbors findKNN
.setup_knn_data <- function(x, subset.row, d, transposed, k, BNPARAM, BSPARAM, BPPARAM) {
    ncells <- ncol(x)
    if (!is.null(subset.row)) {
        x <- x[.subset_to_index(subset.row, x, byrow=TRUE),,drop=FALSE]
    }
    
    if (!transposed) {
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
setGeneric("buildSNNGraph", function(x, ...) standardGeneric("buildSNNGraph"))

#' @export
setMethod("buildSNNGraph", "ANY", .buildSNNGraph)

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim 
#' @export
setMethod("buildSNNGraph", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE, use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        out <- .buildSNNGraph(reducedDim(x, use.dimred), d=NA, transposed=TRUE, ..., subset.row=NULL)
    } else {
        subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)
        out <- .buildSNNGraph(assay(x, i=assay.type), transposed=FALSE, ..., subset.row=subset.row)
    }
    return(out)
})

#' @export
setGeneric("buildKNNGraph", function(x, ...) standardGeneric("buildKNNGraph"))

#' @export
setMethod("buildKNNGraph", "ANY", .buildKNNGraph)

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim 
#' @export
setMethod("buildKNNGraph", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE, use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        out <- .buildKNNGraph(reducedDim(x, use.dimred), d=NA, transposed=TRUE, ..., subset.row=NULL)
    } else {
        subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)
        out <- .buildKNNGraph(assay(x, i=assay.type), transposed=FALSE, ..., subset.row=subset.row)
    }
    return(out)
})
