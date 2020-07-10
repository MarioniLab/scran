#' Compute the cluster-wise Rand indices
#'
#' Breaks down the Rand index calculation to report values for each cluster and pair of clusters.
#'
#' @param ref A character vector or factor containing one set of groupings, considered to be the reference.
#' @param alt A character vector or factor containing another set of groupings, to be compared to \code{alt}.
#' @param mode String indicating whether to return the ratio, the number of pairs or the Rand index.
#' @param adjusted Logical scalar indicating whether the adjusted Rand index should be returned.
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
#' If \code{adjusted=TRUE}, we adjust all counts by subtracting their expected values under a model of random permutations.
#' This accounts for differences in the number and sizes of clusters within and between \code{ref} and \code{alt},
#' in a manner that mimics the calculation of adjusted Rand index (ARI).
#' We subtract expectations on a per-cluster or per-cluster-pair basis for \eqn{a} and \eqn{b}, respectively;
#' we also redefine the \dQuote{total} number of cell pairs for each cluster or cluster pair based on the denominator of the ARI. 
#' 
#' @return
#' If \code{mode="ratio"}, a square numeric matrix is returned with number of rows equal to the number of unique levels in \code{ref}.
#' Each diagonal entry is the ratio of the per-cluster \eqn{a} to the total number of pairs of cells in that cluster.
#' Each off-diagonal entry is the ratio of the per-cluster-pair \eqn{b} to the total number of pairs of cells for that pair of clusters.
#' Lower-triangular entries are set to \code{NA}.
#' If \code{adjusted=TRUE}, counts and totals are both adjusted prior to computing the ratio.
#'
#' If \code{mode="pairs"}, a list is returned containing \code{correct} and \code{total},
#' both of which are square numeric matrices of the same arrangement as described above.
#' However, \code{correct} contains the actual numbers \eqn{a} (diagonal) and \eqn{b} (off-diagonal) rather than the ratios,
#' while \code{total} contains the total number of cell pairs in each cluster or pair of clusters.
#' If \code{adjusted=TRUE}, both matrices are adjusted by subtracting the random expectations from the counts.
#'
#' If \code{mode="index"}, a numeric scalar is returned containing the Rand index (or ARI, if \code{adjusted=TRUE}).
#' 
#' @author Aaron Lun
#' @examples
#' library(scuttle)
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
#' \code{\link{clusterModularity}}, which applies the same breakdown to the cluster modularity.
#'
#' @export
clusterRand <- function(ref, alt, mode=c("ratio", "pairs", "index"), adjusted=TRUE) {
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

    if (adjusted) {
        refn <- table(ref)
        altn <- table(alt)

        # Computing the probability that two randomly chosen cells end up with
        # the same 'alt' label. This is used to compute the expected value for
        # each entry of 'correct'.
        same.alt.n <- sum(choose(altn, 2))
        total.n <- choose(length(alt), 2)
        same.alt.p <- same.alt.n/total.n

        # The LHS of the ARI's denominator is an average of the maximum number
        # of same-cluster pairs in 'ref' and 'alt'. We distribute this number
        # to each cluster in 'ref' according to its number of same-cluster cell
        # pairs. We generalize this for cluster pairs, even though this is not
        # technically used for the calculation of the ARI. 
        same.ref.n <- sum(choose(refn, 2))
        same.alt.mult <- same.alt.n/same.ref.n
        diff.alt.mult <- (total.n - same.alt.n)/(total.n - same.ref.n)

        for (i in seq_along(all.lev)) {
            leftn <- refn[all.lev[i]]
            for (j in seq_len(i)) {
                rightn <- refn[all.lev[j]]
                ref.total <- total[j, i]

                if (i!=j) {
                    expected <- leftn * rightn * (1  - same.alt.p)
                    alt.total <- diff.alt.mult * ref.total
                } else {
                    expected <- choose(leftn, 2) * same.alt.p
                    alt.total <- same.alt.mult * ref.total
                }

                correct[j,i] <- correct[j,i] - expected
                total[j,i] <- 0.5 * (ref.total + alt.total) - expected
            }
        }
    }

    mode <- match.arg(mode)
    if (mode=="ratio") {
        correct/total
    } else if (mode=="index") {
        if (!adjusted) {
            sum(correct, na.rm=TRUE)/sum(total, na.rm=TRUE)
        } else {
            sum(diag(correct))/sum(diag(total))
        }
    } else {
        list(correct=correct, total=total)
    }
}
