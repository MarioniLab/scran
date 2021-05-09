#' Compute the top rank
#'
#' Compute the top rank in a matrix of statistics, usually effect sizes from differential comparisons.
#'
#' @param x A matrix of statistics from multiple differential comparisons (columns) and genes (rows).
#' @param ties.method String specifying how ties should be handled.
#' @param decreasing Logical scalar indicating whether to obtain ranks for decreasing magnitude of values in \code{x}.
#'
#' @details
#' For each gene, the top rank is defined by ranking values within each column of \code{x}, and then taking the smallest rank across rows.
#' This represents the highest rank that each gene achieves in any comparison,
#' assuming that the columns of \code{x} represent statistics from a single comparison.
#'
#' The top rank allows us to easily filter a set of genes to obtain the top DE genes from each comparison.
#' In the context of marker detection with pairwise comparisons between groups of cells, 
#' taking all rows with ranks values less than or equal to T yields a marker set containing the top T genes (ranked by significance) from each pairwise comparison.
#' This guarantees the inclusion of genes that can distinguish between any two groups.
#' 
#' The set of genes with smallest ranks <= 1 will contain the top gene from each pairwise comparison to every other cluster.
#' If T is instead, say, 5, the set will consist of the \emph{union} of the top 5 genes from each pairwise comparison.
#' Multiple genes can have the same smallest rank as different genes may have the same rank across different pairwise comparisons.
#' Conversely, the marker set may be smaller than the product of T and the number of other clusters, as the same gene may be shared across different comparisons.
#'
#' @return A numeric vector containing the smallest (i.e., top) rank for each gene across all comparisons.
#'
#' @seealso
#' \code{\link{scoreMarkers}}, where this function is used to compute one of the effect size summaries.
#' 
#' @examples
#' # Get smallest rank by log-FC:
#' lfcs <- matrix(rnorm(100), ncol=5)
#' computeTopRank(lfcs)
#'
#' # Get smallest rank by p-value: 
#' pvals <- matrix(runif(100), ncol=5)
#' computeTopRank(pvals, decreasing=FALSE)
#' 
#' @export
#' @importFrom DelayedMatrixStats colMins colRanks
computeTopRank <- function(x, ties.method="min", decreasing=TRUE) {
    x <- as.matrix(x)
    if (decreasing) x <- -x
    colMins(colRanks(x, ties.method=ties.method), na.rm=TRUE) 
}
