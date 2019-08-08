#' Compute the cluster-wise modularity
#' 
#' Calculate the modularity of each cluster from a graph, based on a null model of random connections between nodes.
#' 
#' @param graph A \link{graph} object from \pkg{igraph}, usually where each node represents a cell.
#' @param clusters A factor specifying the cluster identity for each node.
#' @param get.weights A logical scalar indicating whether the observed and expected edge weights should be returned, rather than the modularity.
#' @param get.values Deprecated, same as \code{get.weights}.
#' 
#' @return
#' By default, an upper triangular numeric matrix of order equal to the number of clusters is returned.
#' Each entry corresponds to a pair of clusters and is proportional to the difference between the observed and expected edge weights between those clusters.
#' 
#' If \code{as.ratio=TRUE}, an upper triangular numeric matrix is again returned.
#' Here, each entry is equal to the log-ratio between the observed and expected edge weights.
#' 
#' If \code{get.weights=TRUE}, a list is returned containing two upper triangular numeric matrices. 
#' The \code{observed} matrix contains the observed sum of edge weights between and within clusters,
#' while the \code{expected} matrix contains the expected sum of edge weights under the random model.
#' 
#' @details
#' This function computes a modularity score in the same manner as that from \code{\link{modularity}}.
#' The modularity is defined as the (scaled) difference between the observed and expected number of edges between nodes in the same cluster.
#' The expected number of edges is defined by a null model where edges are randomly distributed among nodes.
#' The same logic applies for weighted graphs, replacing the number of edges with the summed weight of edges.
#' 
#' Whereas \code{\link{modularity}} returns a modularity score for the entire graph, \code{clusterModularity} provides scores for the individual clusters.
#' The sum of the diagonal elements of the output matrix should be equal to the output of \code{\link{modularity}} 
#' (after supplying weights to the latter, if necessary).
#' A well-separated cluster should have mostly intra-cluster edges and a high modularity score on the corresponding diagonal entry,
#' while two closely related clusters that are weakly separated will have many inter-cluster edges and a high off-diagonal score.
#'
#' In practice, the modularity may not the most effective metric for evaluating cluster separatedness.
#' This is because the modularity is proportional to the number of cells, so larger clusters will naturally have a large score regardless of separation.
#' An alternative approach is to set \code{as.ratio=TRUE}, which returns the (log-)ratio of the observed to expected weights for each entry of the matrix.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{buildSNNGraph}}, for one method to construct \code{graph}.
#'
#' \code{\link{modularity}}, for the calculation of the entire graph modularity.
#' 
#' @examples
#' example(buildSNNGraph) # using the mocked-up graph in this example.
#' 
#' # Examining the modularity values directly.
#' out <- clusterModularity(g, clusters)
#' image(out)
#' 
#' # Alternatively, compare the ratio of observed:expected.
#' out <- clusterModularity(g, clusters, get.weights=TRUE)
#' log.ratio <- log2(out$observed/out$expected + 1)
#' image(log.ratio)
#' 
#' @export
clusterModularity <- function(graph, clusters, get.values=FALSE) 
# Computes the cluster-wise modularity scores, for use in 
# assessing the quality of graph-based clustering.
#
# written by Aaron Lun
# created 2 January 2018
{
    by.clust <- split(seq_along(clusters), clusters)
    uclust <- names(by.clust)
    nclust <- length(uclust)
    mod.mat <- matrix(0, nclust, nclust)
    rownames(mod.mat) <- colnames(mod.mat) <- uclust

    # Calculating the observed weight within/between clusters.
    grmat <- graph[]
    for (x in seq_along(by.clust)) {
        for (y in seq_len(x)) { 
            current <- by.clust[[x]] 
            other <- by.clust[[y]]
            mod.mat[y,x] <- mod.mat[x,y] <- sum(grmat[current,other])
        }
    }

    # Calcuating the expected weight if they were randomly distributed. 
    total.weight <- sum(mod.mat)
    clust.prop <- colSums(mod.mat)/total.weight
    expected.mat <- tcrossprod(clust.prop) * total.weight

    if (get.values){
        return(list(observed=mod.mat, expected=expected.mat))
    } else {
        return(1/total.weight * (mod.mat - expected.mat))
    }
}

