#' Compute coassignment probabilities
#'
#' Compute coassignment probabilities for each label in a reference grouping when compared to an alternative grouping of samples.
#'
#' @param ref A character vector or factor containing one set of groupings, considered to be the reference.
#' @param alt A character vector or factor containing another set of groupings, to be compared to \code{alt}.
#' @param summarize Logical scalar indicating whether the output matrix should be converted into a per-label summary.
#' 
#' @return 
#' If \code{summarize=FALSE}, a numeric matrix is returned with upper triangular entries filled with the coassignment probabilities for each pair of labels in \code{ref}.
#'
#' Otherwise, a \linkS4class{DataFrame} is returned with one row per label in \code{ref} containing the \code{self} and \code{other} coassignment probabilities.
#'
#' @details
#' The coassignment probability for each pair of labels in \code{ref} is the probability that a randomly chosen cell from each of the two reference labels will have the same label in \code{alt}.
#' High coassignment probabilities indicate that a particular pair of labels in \code{ref} are frequently assigned to the same label in \code{alt}, which has some implications for cluster stability.
#'
#' When \code{summarize=TRUE}, we summarize the matrix of coassignment probabilities into a set of per-label values.
#' The \dQuote{self} coassignment probability is simply the diagonal entry of the matrix, i.e., the probability that two cells from the same label in \code{ref} also have the same label in \code{alt}.
#' The \dQuote{other} coassignment probability is the maximum probability across all pairs involving that label.
#'
#' % One might consider instead reporting the 'other' probability as the probability that a randomly chosen cell in the cluster and a randomly chosen cell in any other cluster belong in the same cluster.
#' % However, this results in very small probabilities in all cases, simply because most of the other clusters are well seperated.
#' % Reporting the maximum is more useful as at least you can tell that a cluster is well-separated from _all_ other clusters if it has a low 'other' probability.
#'
#' In general, \code{ref} is well-recapitulated by \code{alt} if the diagonal entries of the matrix is much higher than the sum of the off-diagonal entries.
#' This manifests as higher values for the self probabilities compared to the other probabilities.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{bootstrapCluster}}, to compute coassignment probabilities across bootstrap replicates.
#'
#' \code{\link{clusterRand}}, for another way to compare different clusterings.
#'
#' @examples
#' library(scater)
#' sce <- mockSCE(ncells=200)
#' sce <- logNormCounts(sce)
#' 
#' clust1 <- kmeans(t(logcounts(sce)),3)$cluster
#' clust2 <- kmeans(t(logcounts(sce)),5)$cluster
#'
#' coassignProb(clust1, clust2)
#' coassignProb(clust1, clust2, summarize=TRUE)
#' 
#' @export
coassignProb <- function(ref, alt, summarize=FALSE) {
    tab <- table(ref, alt)
    nref <- nrow(tab)
    output <- matrix(NA_real_, nref, nref,
        dimnames=list(rownames(tab), rownames(tab)))

    for (j1 in seq_len(nref)) {
        spread1 <- tab[j1,]
        spread1 <- spread1/sum(spread1)

        for (j2 in seq_len(j1)) {
            spread2 <- tab[j2,]
            spread2 <- spread2/sum(spread2)

            output[j2,j1] <- sum(spread1 * spread2)
        }
    }

    if (!summarize) {
        output
    } else {
        .summarize_coassign(output)
    }
}

#' @importFrom Matrix forceSymmetric rowSums
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray rowMaxs
.summarize_coassign <- function(mat) {
    flipped <- forceSymmetric(mat)
    self <- diag(flipped)
    diag(flipped) <- 0

    # Define against NA coassignments when empty levels are involved.
    DataFrame(self=self, other=rowMaxs(as.matrix(flipped), na.rm=TRUE), 
        row.names=rownames(flipped))
}
