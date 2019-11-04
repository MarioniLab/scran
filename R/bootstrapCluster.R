#' Assess cluster stability by bootstrapping
#'
#' Generate bootstrap replicates and recluster on them to determine the stability of clusters with respect to sampling noise.
#'
#' @param x A two-dimensional object containing cells in columns.
#' This is usually a numeric matrix of log-expression values, but can also be a \linkS4class{SingleCellExperiment} object.
#' 
#' If \code{transposed=TRUE}, cells are expected to be in the rows, e.g., for precomputed PCs.
#' @param FUN A function that accepts \code{x} and returns a vector of cluster identities.
#' @param clusters A vector of cluster identities obtained by calling \code{FUN(x, ...)}.
#' This is provided as an additional argument in the case that the clusters have already been computed,
#' in which case we can save a single round of computation.
#' @param transposed Logical scalar indicating whether \code{x} is transposed with cells in the rows.
#' @param iterations Integer scalar specifying the number of bootstrap iterations.
#' @param ... Further arguments to pass to \code{FUN} to control the clustering procedure.
#'
#' @return A numeric matrix with upper triangular entries filled with the co-assignment probabilities for each pair of clusters in \code{clusters}.
#'
#' @details
#' Bootstrapping is conventionally used to evaluate the precision of an estimator by applying it to an \emph{in silico}-generated replicate dataset.
#' We can (ab)use this framework to determine the stability of the clusters in the context of a scRNA-seq analysis.
#' We sample cells with replacement from \code{x}, perform clustering with \code{FUN} and compare the new clusters to \code{clusters}.
#'
#' The relevant statistic is the co-assignment probability for each pair of original clusters, i.e., the probability that a randomly chosen cells from each of the two original clusters will be put in the same bootstrap cluster.
#' High co-assignment probabilities indicate that the two original clusters were not stably separated.
#' We might then only trust separation between two clusters if their co-assignment probability was less than some threshold, e.g., 5\%.
#'
#' @section Statistical notes:
#' We use the co-assignment probability as it is more interpretable than, say, the Jaccard index (see the \pkg{fpc} package).
#' It also focuses on the relevant differences between clusters, allowing us to determine which aspects of a clustering are stable.
#' For example, A and B may be well separated but A and C may not be, which is difficult to represent in a single stability measure for A.
#' If our main interest lies in the A/B separation, we do not want to be overly pessimistic about the stability of A, even though it might not be well-separated from all other clusters.
#' 
#' Technically speaking, some mental gymnastics are required to compare the original and bootstrap clusters in this manner.
#' After bootstrapping, the sampled cells represent distinct entities from the original dataset (otherwise it would be difficult to treat them as independent replicates) for which the original clusters do not immediately apply.
#' Instead, we assume that we perform label transfer using a nearest-neighbors approach - which, in this case, is the same as using the original label for each cell, as the nearest neighbor of each resampled cell to the original dataset is itself.
#'
#' Needless to say, bootstrapping will only generate replicates that differ by sampling noise.
#' Real replicates will differ due to composition differences, variability in expression across individuals, etc.
#' Thus, any stability inferences from bootstrapping are likely to be overly optimistic.
#' 
#' @author Aaron Lun
#' @examples
#' library(scater)
#' sce <- mockSCE(ncells=200)
#' bootstrapCluster(sce, FUN=quickCluster, min.size=10)
#'
#' # Defining your own function:
#' sce <- logNormCounts(sce)
#' sce <- runPCA(sce)
#' bootstrapCluster(reducedDim(sce), transposed=TRUE, FUN=function(x) {
#'     kmeans(x, 2)$cluster  
#' })
#'
#' @seealso
#' \code{\link{quickCluster}}, to get a quick and dirty function to use as \code{FUN}.
#' It is often more computationally efficient to define your own function, though.
#' @export
bootstrapCluster <- function(x, FUN, clusters=NULL, transposed=FALSE, iterations=20, ...) {
    if (is.null(clusters)) {
        clusters <- FUN(x, ...)
    }

    cluster.ids <- as.character(sort(unique(clusters)))
    output <- matrix(0, length(cluster.ids), length(cluster.ids))
    output[lower.tri(output, diag=TRUE)] <- NA_real_
    dimnames(output) <- list(cluster.ids, cluster.ids)

    for (i in seq_len(iterations)) {
        if (transposed) {
            chosen <- sample(nrow(x), nrow(x), replace=TRUE)
            resampled <- x[chosen,,drop=FALSE]
        } else {
            chosen <- sample(ncol(x), ncol(x), replace=TRUE)
            resampled <- x[,chosen,drop=FALSE]
        }

        reclusters <- FUN(resampled, ...)
        tab <- table(clusters[chosen], reclusters)

        for (j1 in seq_along(cluster.ids)) {
            spread1 <- tab[cluster.ids[j1],]
            spread1 <- spread1/sum(spread1)

            for (j2 in seq_len(j1-1L)) {
                spread2 <- tab[cluster.ids[j2],]
                spread2 <- spread2/sum(spread2)
                output[j2,j1] <- output[j2,j1] + sum(spread1 * spread2)/iterations
            }
        }
    }

    output
}
