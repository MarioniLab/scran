#' Compute the cluster-wise modularity
#' 
#' Calculate the modularity of each cluster from a graph, based on a null model of random connections between nodes.
#' 
#' @param graph A \link{graph} object from \pkg{igraph}, usually where each node represents a cell.
#' @param clusters Factor specifying the cluster identity for each node.
#' @param get.weights Logical scalar indicating whether the observed and expected edge weights should be returned, rather than the modularity.
#' @param as.ratio Logical scalar indicating whether the log-ratio of observed to expected weights should be returned.
#' 
#' @return
#' By default, an upper triangular numeric matrix of order equal to the number of clusters is returned.
#' Each entry corresponds to a pair of clusters and is proportional to the difference between the observed and expected edge weights between those clusters.
#' 
#' If \code{as.ratio=TRUE}, an upper triangular numeric matrix is again returned.
#' Here, each entry is equal to the ratio between the observed and expected edge weights.
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
#' An alternative approach is to set \code{as.ratio=TRUE}, which returns the ratio of the observed to expected weights for each entry of the matrix.
#' This adjusts for differences in cluster size and improves resolution of differences between clusters.
#'
#' Directed graphs are treated as undirected inputs with \code{mode="each"} in \code{\link{as.undirected}}.
#' In the rare case that self-loops are present, these will also be handled correctly.
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{buildSNNGraph}}, for one method to construct \code{graph}.
#'
#' \code{\link{modularity}}, for the calculation of the entire graph modularity.
#'
#' \code{\link{clusterRand}}, which applies a similar breakdown to the Rand index.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' g <- buildSNNGraph(sce)
#' clusters <- igraph::cluster_walktrap(g)$membership
#' 
#' # Examining the modularity values directly.
#' out <- clusterModularity(g, clusters)
#' out
#' 
#' # Compute the ratio instead, for visualization
#' # (log-transform to improve range of colors).
#' out <- clusterModularity(g, clusters, as.ratio=TRUE)
#' image(log2(out+1))
#'
#' # This can also be used to construct a graph of clusters,
#' # for use in further plotting, a.k.a. graph abstraction.
#' # (Fiddle with the scaling values for a nicer plot.)
#' g2 <- igraph::graph_from_adjacency_matrix(out, mode="upper",
#'     diag=FALSE, weighted=TRUE)
#' plot(g2, edge.width=igraph::E(g2)$weight*10,
#'     vertex.size=sqrt(table(clusters))*10)
#' 
#' # Alternatively, get the edge weights directly:
#' out <- clusterModularity(g, clusters, get.weights=TRUE)
#' out
#'
#' @export
#' @importFrom Matrix diag diag<-
#' @importFrom igraph is.directed
clusterModularity <- function(graph, clusters, get.weights=FALSE, as.ratio=FALSE) {
    by.clust <- split(seq_along(clusters), clusters)
    uclust <- names(by.clust)
    nclust <- length(uclust)

    mod.mat <- matrix(0, nclust, nclust)
    dimnames(mod.mat) <- list(uclust, uclust)
    clust.total <- numeric(nclust)
    names(clust.total) <- uclust

    fullmat <- graph[]

    # Calculating the observed weight within/between clusters.
    for (x in seq_along(by.clust)) {
        current <- by.clust[[x]] 
        for (y in seq_len(x)) { 
            other <- by.clust[[y]]
            grmat <- fullmat[current,other,drop=FALSE]
            grsum <- sum(grmat)
                
            if (x==y) {
                old.diag <- sum(diag(grmat))
                diag(grmat) <- 0

                self.sum <- sum(grmat)
                if (!is.directed(graph)) {
                    # Only count undirected edges within a cluster once. 
                    self.sum <- self.sum/2
                } else {
                    # Need to count directed edges between different nodes twice.
                    grsum <- grsum + self.sum
                }

                # Self-edges contribute twice to total node weight,
                # according to igraph::modularity.
                grsum <- grsum + old.diag

                mod.mat[x,y] <- self.sum + old.diag 
            } else {
                if (is.directed(graph)) {
                    grsum <- grsum + sum(graph[other,current])
                }
                mod.mat[y,x] <- grsum
            }

            # If x==y, this is equivalent to adding edge weights for each node twice.
            # THIS IS DELIBERATE; the total per-node weight is that for all edges 
            # involving each node, so edges between nodes in the same cluster are 
            # allowed to be double-counted here when aggregating node weights per cluster.
            clust.total[x] <- clust.total[x] + grsum
            if (x!=y) {
                clust.total[y] <- clust.total[y] + grsum
            }
        }
    }

    # Calcuating the expected weight if they were randomly distributed. 
    # Note some effort is involved in 'folding' it into an upper triangular.
    total.weight <- sum(mod.mat)
    clust.prop <- clust.total/sum(clust.total)
    expected.mat <- tcrossprod(clust.prop) * total.weight

    expected.mat[lower.tri(expected.mat)] <- 0
    old.diag <- diag(expected.mat)
    expected.mat <- expected.mat * 2
    diag(expected.mat) <- old.diag

    if (get.weights) {
        list(observed=mod.mat, expected=expected.mat)
    } else if (as.ratio) {
        mod.mat/expected.mat
    } else {
        1/total.weight * (mod.mat - expected.mat)
    }
}
