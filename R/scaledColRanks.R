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
    return(rkout)
}
