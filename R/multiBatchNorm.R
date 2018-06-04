#' @export
#' @importFrom scater calcAverage
#' @importFrom BiocGenerics normalize sizeFactors sizeFactors<-
#' @importFrom stats median
multiBatchNorm <- function(..., assay.type="counts", norm.args=list(), subset.row=NULL, min.mean=1)
# Performs multi-batch normalization, adjusting for differences 
# in scale between SCE objects supplied in '...'.
# 
# written by Aaron Lun
# created 4 June 2018
{
    batches <- list(...)
    nbatches <- length(batches)
    collected.ave <- vector("list", nbatches)
    for (idx in seq_along(batches)) {
        collected.ave[[idx]] <- calcAverage(batches[[idx]], subset_row=subset.row, exprs_values=assay.type, use_size_factors=FALSE)
    }

    # Computing the median ratios (a la DESeq normalization).
    collected.ratios <- matrix(1, nbatches, nbatches)
    for (first in seq_len(nbatches-1L)) {
        first.ave <- collected.ave[[first]]
        first.sum <- sum(first.ave)
        for (second in first + seq_len(nbatches - first)) {
            second.ave <- collected.ave[[second]]
            second.sum <- sum(second.ave)

            # Mimic calcAverage(cbind(first.ave, second.ave)).
            grand.mean <- (first.ave/first.sum + second.ave/second.sum)/2 * (first.sum + second.sum)/2
            keep <- grand.mean >= min.mean

            curratio <- median(second.ave[keep]/first.ave[keep])
            collected.ratios[first,second] <- curratio
            collected.ratios[second,first] <- 1/curratio
        }
    }

    # Finding the smallest ratio, and using that as the reference.
    smallest <- which.min(apply(collected.ratios, 2, min))
    
    for (idx in seq_along(batches)) {
        current <- batches[[idx]]
        cursf <- sizeFactors(current)
        if (is.null(cursf)) {
            cursf <- librarySizeFactors(current)
            warning("using library sizes as size factors")
        }
    
        # Adjusting the size factors for the new reference.
        cursf <- cursf/mean(cursf) / collected.ratios[idx,smallest]
        sizeFactors(current) <- cursf
        batches[[idx]] <- do.call(normalize, c(list(current, exprs_values=assay.type, centre_size_factors=FALSE), norm.args))
    }

    return(batches)
}
