#' @export
#' @importFrom scater librarySizeFactors
#' @importFrom SingleCellExperiment sizeFactorNames
#' @importFrom BiocGenerics sizeFactors sizeFactors<- normalize
#' @importFrom stats ave
multiBlockNorm <- function(x, block, ...) 
# Adjusts size factors so that the spike-in size factors have
# the same mean as the reference size factors _within each batch_.
# Ensures that the abundances are comparable in multiBlockVar.
#
# written by Aaron Lun
# created 4 June 2018
{
    ref <- sizeFactors(x)
    if (is.null(ref)) {
        ref <- librarySizeFactors(x)
        warning("using library sizes as size factors")
    } else {
        # Centering just in case.
        ref <- ref / mean(ref)
        sizeFactors(x) <- ref
    }
    ref.mean <- ave(ref, block, FUN=mean)

    for (sf in sizeFactorNames(x)) {
        current <- sizeFactors(x, sf)
        cur.mean <- ave(current, block, FUN=mean)
        sizeFactors(x, sf) <- current * ref.mean / cur.mean
    }
   
    x <- normalize(x, ..., centre_size_factors=FALSE)
    return(x)
}
