#' Assess cluster stability by bootstrapping
#'
#' Generate bootstrap replicates and recluster on them to determine the stability of clusters with respect to sampling noise.
#'
#' @param x A two-dimensional object containing cells in columns.
#' This is usually a numeric matrix of log-expression values, 
#' but can also be a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object.
#' 
#' If \code{transposed=TRUE}, a matrix is expected where cells are in the rows, e.g., for precomputed PCs.
#' @param FUN A function that accepts a value of the same type as \code{x} and returns a vector or factor of cluster identities.
#' @param clusters A vector or factor of cluster identities equivalent to that obtained by calling \code{FUN(x, ...)}.
#' This is provided as an additional argument in the case that the clusters have already been computed,
#' in which case we can save a single round of computation.
#' @param transposed Logical scalar indicating whether \code{x} is transposed with cells in the rows.
#' @param iterations A positive integer scalar specifying the number of bootstrap iterations.
#' @param average String specifying the method to use to average across bootstrap iterations.
#' @param ... Further arguments to pass to \code{FUN} to control the clustering procedure.
#' @param compare A function that accepts the original clustering and the bootstrapped clustering,
#' and returns a numeric vector or matrix containing some measure of similarity between them - see Details.
#' @param mode,adjusted Further arguments to pass to \code{\link{clusterRand}} when \code{compare=NULL}.
#'
#' @return 
#' If \code{compare=NULL} and \code{mode="ratio"}, a numeric matrix is returned with upper triangular entries filled with the ratio of the adjusted cell pair counts (see \code{?\link{clusterRand}}) for each pair of clusters in \code{clusters}.
#' Each ratio is averaged across bootstrap iterations.
#' 
#' If \code{compare=NULL} and \code{mode="index"}, a numeric scalar containing the average ARI between \code{clusters} and the bootstrap replicates across iterations is returned.
#'
#' If \code{compare} is provided, a numeric array of the same type as the output of \code{compare} is returned,
#' containing the average statistic(s) across bootstrap replicates.
#'
#' @details
#' Bootstrapping is conventionally used to evaluate the precision of an estimator by applying it to an \emph{in silico}-generated replicate dataset.
#' We can (ab)use this framework to determine the stability of the clusters in the context of a scRNA-seq analysis.
#' We sample cells with replacement from \code{x}, perform clustering with \code{FUN} and compare the new clusters to \code{clusters}.
#'
#' For comparing clusters, we compute the ratio matrix from \code{\link{clusterRand}} and average its values across bootstrap iterations.
#' High on-diagonal values indicate that the corresponding cluster remains coherent in the bootstrap replicates,
#' while high off-diagonal values indicate that the corresponding pair of clusters are still separated in the replicates.
#' If a single value is necessary, we can instead average the adjusted Rand indices across iterations with \code{mode="ratio"}.
#' 
#' We use the ratio matrix by default as it is more interpretable than a single value like the ARI or the Jaccard index (see the \pkg{fpc} package).
#' It focuses on the relevant differences between clusters, allowing us to determine which aspects of a clustering are stable.
#' For example, A and B may be well separated but A and C may not be, which is difficult to represent in a single stability measure for A.
#' If our main interest lies in the A/B separation, we do not want to be overly pessimistic about the stability of A, even though it might not be well-separated from all other clusters.
#'
#' @section Using another comparison function:
#' We can use a different method for comparing clusterings by setting \code{compare}.
#' This is expected to be a function that takes two arguments - 
#' the original clustering first, and the bootstrapped clustering second - 
#' and returns some kind of numeric scalar, vector or matrix containing 
#' statistics for the similarity or difference between the original and bootstrapped clustering.
#' These statistics are then averaged across all bootstrap iterations.
#' 
#' Any numeric output of \code{compare} is acceptable as long as the dimensions are only dependent on the \emph{levels} of the original clustering - including levels that have no cells, due to resampling! - and thus do not change across bootstrap iterations.
#' One example of a compatible function is \code{\link{clusterRand}}, 
#' which provides a cluster-wise breakdown of the Rand index.
#'
#' @section Statistical note on bootstrap comparisons:
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
#' library(scuttle)
#' sce <- mockSCE(ncells=200)
#'
#' # Using 'quickCluster' as one potential 'FUN':
#' bootstrapCluster(sce, FUN=quickCluster, min.size=10)
#'
#' # Defining your own function:
#' sce <- logNormCounts(sce)
#' sce <- scater::runPCA(sce)
#' kFUN <- function(x) kmeans(x, 2)$cluster  
#' bootstrapCluster(reducedDim(sce), transposed=TRUE, FUN=kFUN)
#'
#' # Using an alternative comparison, in this case the Rand index:
#' bootstrapCluster(reducedDim(sce), transposed=TRUE, FUN=kFUN, 
#'     compare=clusterRand)
#'
#' @seealso
#' \code{\link{quickCluster}}, to get a quick and dirty function to use as \code{FUN}.
#' It is often more computationally efficient to define your own function, though.
#'
#' \code{\link{clusterRand}}, for the calculation of the ARI.
#' 
#' @export
#' @importFrom DelayedMatrixStats rowMedians
bootstrapCluster <- function(x, FUN, clusters=NULL, transposed=FALSE, iterations=20, 
    average=c("median", "mean"), ..., compare=NULL, mode="ratio", adjusted=TRUE)
{
    if (is.null(clusters)) {
        clusters <- FUN(x, ...)
    }
    clusters <- as.factor(clusters)

    if (iterations <= 0L) {
        stop("'iterations' must be a positive integer")
    }

    collated <- vector("list", iterations)
    if (is.null(compare)) {
        compare <- function(...) clusterRand(..., mode=mode, adjusted=adjusted)
    }

    for (i in seq_len(iterations)) {
        if (transposed) {
            chosen <- sample(nrow(x), nrow(x), replace=TRUE)
            resampled <- x[chosen,,drop=FALSE]
        } else {
            chosen <- sample(ncol(x), ncol(x), replace=TRUE)
            resampled <- x[,chosen,drop=FALSE]
        }

        reclusters <- FUN(resampled, ...)
        collated[[i]] <- compare(clusters[chosen], reclusters)
    }

    # A robust way of computing the average that handles NAs.
    if (length(unique(lapply(collated, dim))) > 1L) { 
        stop("'compare' output should have constant dimension")
    }

    as.mat <- do.call(cbind, lapply(collated, as.numeric))

    average <- match.arg(average)
    if (average=="mean") {
        averaged <- rowMeans(as.mat, na.rm=TRUE)
    } else {
        averaged <- rowMedians(as.mat, na.rm=TRUE)
    }

    dim(averaged) <- dim(collated[[1]])
    dimnames(averaged) <- dimnames(collated[[1]])
    averaged 
}
