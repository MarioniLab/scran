#' @importFrom igraph make_graph E simplify "E<-"
#' @importFrom BiocParallel SerialParam
.buildSNNGraph <- function(x, k=10, d=50, transposed=FALSE, pc.approx=FALSE,
                           rand.seed=NA, irlba.args=list(), 
                           subset.row=NULL, BPPARAM=SerialParam()) 
# Builds a shared nearest-neighbor graph, where edges are present between each 
# cell and any other cell with which it shares at least one neighbour. Each edges 
# is weighted based on the ranks of the shared nearest neighbours of the two cells, 
# as described in the SNN-Cliq paper.
#
# written by Aaron Lun
# created 3 April 2017
# last modified 16 November 2017    
{ 
    nn.out <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed,
        pc.approx=pc.approx, rand.seed=rand.seed, irlba.args=irlba.args, 
        k=k, BPPARAM=BPPARAM) 

    # Building the SNN graph.
    g.out <- .Call(cxx_build_snn, nn.out$index)
    edges <- g.out[[1]] 
    weights <- g.out[[2]]

    g <- make_graph(edges, directed=FALSE)
    E(g)$weight <- weights
    g <- simplify(g, edge.attr.comb="first") # symmetric, so doesn't really matter.
    return(g)
}

#' @importFrom igraph make_graph simplify
#' @importFrom BiocParallel SerialParam
.buildKNNGraph <- function(x, k=10, d=50, directed=FALSE, transposed=FALSE, pc.approx=FALSE,
                           rand.seed=NA, irlba.args=list(), 
                           subset.row=NULL, BPPARAM=SerialParam()) 
# Builds a k-nearest-neighbour graph, where edges are present between each
# cell and its 'k' nearest neighbours. Undirected unless specified otherwise.
#
# written by Aaron Lun, Jonathan Griffiths
# created 16 November 2017
{ 
    nn.out <- .setup_knn_data(x=x, subset.row=subset.row, d=d, transposed=transposed,
        pc.approx=pc.approx, rand.seed=rand.seed, irlba.args=irlba.args,
        k=k, BPPARAM=BPPARAM) 

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
#' @importFrom kmknn findKNN
.setup_knn_data <- function(x, subset.row, d, transposed, pc.approx, rand.seed, irlba.args, k, BPPARAM) {
    ncells <- ncol(x)
    if (!is.null(subset.row)) {
        x <- x[.subset_to_index(subset.row, x, byrow=TRUE),,drop=FALSE]
    }
    
    if (!transposed) {
        x <- t(x)
    } 
    
    # Reducing dimensions, if 'd' is less than the number of genes.
    if (!is.na(d) && d < ncol(x)) {
        if (pc.approx) {
            if (!is.na(rand.seed)) {
                .Deprecated(msg="'rand.seed=' is deprecated.\nUse 'set.seed' externally instead.")
                set.seed(rand.seed)
            }
            pc <- do.call(irlba::prcomp_irlba, c(list(x=x, n=d, scale.=FALSE, center=TRUE, retx=TRUE), irlba.args))
        } else {
            pc <- prcomp(x, rank.=d, scale.=FALSE, center=TRUE)
        }
        x <- pc$x
    }
   
    # Finding the KNNs. 
    findKNN(x, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
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
setMethod("buildSNNGraph", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE, use.dimred=NULL) {
              
    subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)
    if (!is.null(use.dimred)) {
        out <- .buildSNNGraph(reducedDim(x, use.dimred), d=NA, transposed=TRUE, ..., subset.row=NULL)
    } else {
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
setMethod("buildKNNGraph", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE, use.dimred=NULL) {
              
    subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)
    if (!is.null(use.dimred)) {
        out <- .buildKNNGraph(reducedDim(x, use.dimred), d=NA, transposed=TRUE, ..., subset.row=NULL)
    } else {
        out <- .buildKNNGraph(assay(x, i=assay.type), transposed=FALSE, ..., subset.row=subset.row)
    }
    return(out)
})
