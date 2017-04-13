.buildSNNGraph <- function(x, k=10, d=50, subset.row=NULL) 
# Builds a shared nearest-neighbor graph, where edges are present between each 
# cell and its 'k' nearest neighbours. Edges are weighted based on the ranks of 
# the shared nearest neighbours of the two cells, as described in the SNN-Cliq paper.
#
# written by Aaron Lun
# created 3 April 2017
# last modified 13 April 2017    
{ 
    ncells <- ncol(x)
    if (!is.null(subset.row)) {
        x <- x[.subset_to_index(subset.row, x, byrow=TRUE),,drop=FALSE]
    }
    
    # Reducing dimensions.
    x <- t(x)
    if (!is.na(d) && d < ncells) {
        pc <- prcomp(x)
        x <- pc$x[,seq_len(d),drop=FALSE]
    }

    # Getting the kNNs.
    nn.out <- get.knn(x, k=k, algorithm="cover_tree") 

    g.out <- .Call(cxx_build_snn, nn.out$nn.index)
    if (is.character(g.out)) { stop(g.out) }
    edges <- g.out[[1]] 
    weights <- g.out[[2]]

    g <- make_graph(edges, directed=FALSE)
    E(g)$weight <- weights
    g <- simplify(g, edge.attr.comb="first") # symmetric, so doesn't really matter.
    return(g)
}

setGeneric("buildSNNGraph", function(x, ...) standardGeneric("buildSNNGraph"))

setMethod("buildSNNGraph", "matrix", .buildSNNGraph)

setMethod("buildSNNGraph", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { 
        subset.row <- .spikeSubset(x, get.spikes)
    }
    .buildSNNGraph(assayDataElement(x, assay), ..., subset.row=subset.row)
})
