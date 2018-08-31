#' @importFrom scater calcAverage
.scaled_col_ranks <- function(x, subset.row=NULL, min.mean=NULL, transposed=FALSE, withDimnames=TRUE)
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

    rkout <- .Call(cxx_get_scaled_ranks, x, subset.row-1L, transposed)
    if (withDimnames) {
        dimnames(rkout) <- dimnames(x)
    }
    return(rkout)
}

#' @export
setGeneric("scaledColRanks", function(x, ...) standardGeneric("scaledColRanks"))

#' @export
setMethod("scaledColRanks", "ANY", .scaled_col_ranks)

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("scaledColRanks", "SingleCellExperiment", function(x, subset.row=NULL, ..., assay.type="counts", get.spikes=FALSE) { 
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)          
    .scaled_col_ranks(assay(x, i=assay.type), subset.row=subset.row, ...)
})
