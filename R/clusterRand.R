#' Compute the cluster-wise Rand indices
#'
#' Breaks down the Rand index calculation to report values for each cluster and pair of clusters.
#'
#' @inheritParams coassignProb
#' @param mode String indicating whether to return the ratio, the number of pairs or the Rand index.
#'
#' @details
#' Recall that the Rand index calculation consists of four numbers:
#' \describe{
#' \item{\eqn{a}}{The number of pairs of cells in the same cluster in \code{ref} and the same cluster in \code{alt}.}
#' \item{\eqn{b}}{The number of pairs of cells in different clusters in \code{ref} and different clusters in \code{alt}.}
#' \item{\eqn{c}}{The number of pairs of cells in the same cluster in \code{ref} and different clusters in \code{alt}.}
#' \item{\eqn{d}}{The number of pairs of cells in different clusters in \code{ref} but the same cluster in \code{alt}.}
#' }
#' The Rand index is then computed as \eqn{a + b} divided by \eqn{a + b + c + d}, i.e., the total number of pairs.
#'
#' We can break these numbers down into values for each cluster or pair of clusters in \code{ref}.
#' For each cluster, we compute its value of \eqn{a}, 
#' i.e., the number of pairs of cells in \emph{that} cluster that are also in the same cluster in \code{alt}.
#' Similarly, for each pair of clusters in \code{ref}, we compute its value of \eqn{b},
#' i.e., the number of pairs of cells that have one cell in each of those clusters 
#' and also belong in different clusters in \code{alt}.
#'
#' This process provides more information about the specific similarities or differences between \code{ref} and \code{alt},
#' rather than coalescing all the values into a single statistic. 
#' For example, it is now possible to see which specific clusters from \code{ref} are not reproducible in \code{alt},
#' or which specific partitions between pairs of clusters are not reproducible.
#' In the default output, such events can be diagnosed by looking for low entries in the ratio matrix;
#' on the other hand, values close to 1 indicate that \code{ref} is almost perfectly recapitulated by \code{alt}.
#' 
#' @return
#' By default, a square numeric matrix with number of rows equal to the number of unique levels in \code{ref}.
#' Each diagonal entry represents the ratio of the per-cluster \eqn{a} to the total number of pairs of cells in that cluster.
#' Each off-diagonal entry represents the ratio of the per-cluster-pair \eqn{b} to the total number of pairs of cells for that pair of clusters.
#' Lower-triangular entries are set to \code{NA}.
#'
#' If \code{mode="pairs"}, a list is returned containing \code{correct} and \code{total},
#' both of which are square numeric matrices of the same arrangement as described above.
#' However, \code{correct} contains the actual numbers \eqn{a} (diagonal) and \eqn{b} (off-diagonal) rather than the ratios,
#' while \code{total} contains the total number of cell pairs in each cluster or pair of clusters.
#'
#' If \code{mode="index"}, a numeric scalar is returned containing the Rand index.
#' 
#' @author Aaron Lun
#' @examples
#' library(scater)
#' sce <- mockSCE(ncells=200)
#' sce <- logNormCounts(sce)
#' 
#' clust1 <- kmeans(t(logcounts(sce)),3)$cluster
#' clust2 <- kmeans(t(logcounts(sce)),5)$cluster
#'
#' ratio <- clusterRand(clust1, clust2)
#' ratio
#'
#' # Getting the raw counts:
#' clusterRand(clust1, clust2, mode="pairs")
#' 
#' # Computing the original Rand index.
#' clusterRand(clust1, clust2, mode="index")
#'
#' @seealso
#' \code{\link{coassignProb}}, for another way of comparing two clusterings.
#'
#' \code{\link{clusterModularity}}, which applies the same breakdown to the cluster modularity.
#'
#' @export
clusterRand <- function(ref, alt, mode=c("ratio", "pairs", "index")) {
    ref <- as.factor(ref)
    alt <- as.factor(alt)
    all.lev <- levels(ref)

    correct <- total <- matrix(NA_real_, length(all.lev), length(all.lev), 
        dimnames=list(all.lev, all.lev))

    for (i in seq_along(all.lev)) {
        left <- alt[ref==all.lev[i]]
        tab.left <- table(left)

        for (j in seq_len(i)) {
            right <- alt[ref==all.lev[j]]
            tab.right <- table(right)

            if (i!=j) {
                total.pairs <- sum(tab.left) * sum(tab.right) 
                correct.pairs <- total.pairs - sum(tab.left * tab.right)
            } else {
                correct.pairs <- sum(choose(tab.left, 2))
                total.pairs <- choose(sum(tab.left), 2)
            }

            correct[j,i] <- correct.pairs
            total[j,i] <- total.pairs
        }
    }

    mode <- match.arg(mode)
    if (mode=="ratio") {
        correct/total
    } else if (mode=="index") {
        sum(correct, na.rm=TRUE)/sum(total, na.rm=TRUE)
    } else {
        list(correct=correct, total=total)
    }
}
