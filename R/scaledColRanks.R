#' Compute scaled column ranks
#'
#' Compute scaled column ranks from each cell's expression profile for distance calculations based on rank correlations.
#' 
#' @param x A numeric matrix-like object containing cells in columns and features in the rows.
#' @param subset.row A logical, integer or character scalar indicating the rows of \code{x} to use, see \code{?"\link{scran-gene-selection}"}.
#' @param min.mean A numeric scalar specifying the filter to be applied on the average normalized count for each feature prior to computing ranks.
#' Disabled by setting to \code{NULL}.
#' @param transposed A logical scalar specifying whether the output should be transposed.
#' @param as.sparse A logical scalar indicating whether the output should be sparse.
#' @param withDimnames A logical scalar specifying whether the output should contain the dimnames of \code{x}.
#' 
#' @return
#' A matrix of the same dimensions as \code{x}, where each column contains the centred and scaled ranks of the expression values for each cell.
#' If \code{transposed=TRUE}, this matrix is transposed so that rows correspond to cells.
#' If \code{as.sparse}, the columns are not centered to preserve sparsity.
#' 
#' @details
#' Euclidean distances computed based on the output rank matrix are equivalent to distances computed from Spearman's rank correlation.
#' This can be used in clustering, nearest-neighbour searches, etc. as a robust alternative to Euclidean distances computed directly from \code{x}. 
#' 
#' If \code{as.sparse=TRUE}, the most common average rank is set to zero in the output.
#' This can be useful for highly sparse input data where zeroes have the same rank and are themselves returned as zeroes.
#' Obviously, this means that the ranks are not centred, so this will have to be done manually prior to any downstream distance calculations.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{quickCluster}}, where this function is used.
#' 
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' rout <- scaledColRanks(counts(sce), transposed=TRUE)
#' 
#' # For use in clustering:
#' d <- dist(rout)
#' table(cutree(hclust(d), 4))
#' 
#' g <- buildSNNGraph(rout, transposed=TRUE)
#' table(igraph::cluster_walktrap(g)$membership)
#' 
#' @export
#' @importFrom scater calcAverage
scaledColRanks <- function(x, subset.row=NULL, min.mean=NULL, transposed=FALSE, as.sparse=FALSE, withDimnames=TRUE)
# Obtaining scaled/centred ranks to compute cosine distances.
# Using this instead of colRanks to support dgCMatrix, HDF5Array objects.
# 
# written by Aaron Lun
# created 31 August 2018
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    if (!is.null(min.mean) && all(dim(x)>0L)) {
        further.subset <- calcAverage(x, subset_row=subset.row) >= min.mean
        subset.row <- subset.row[further.subset]
    }

    rkout <- get_scaled_ranks(x, subset.row-1L, transposed, as.sparse)

    if (withDimnames && !is.null(dimnames(x))) {
        dn <- list(rownames(x)[subset.row], colnames(x))
        if (transposed) { 
            dn <- rev(dn)
        }
        dimnames(rkout) <- dn
    }
    rkout
}
