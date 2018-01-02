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

