#' @export
#' @importFrom scater calcAverage
#' @importFrom BiocGenerics normalize sizeFactors sizeFactors<-
#' @importFrom stats median
multiBatchNorm <- function(..., assay.type="counts", norm.args=list(), min.mean=1)
# Performs multi-batch normalization, adjusting for differences 
# in scale between SCE objects supplied in '...'.
# 
# written by Aaron Lun
# created 4 June 2018
{
    batches <- list(...)
    nbatches <- length(batches)
    .check_batch_consistency(batches, byrow=TRUE, ignore.null=TRUE)

    # Computing the median ratios (a la DESeq normalization).
    collected.ave <- lapply(batches, calcAverage, use_size_factors=FALSE)
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
            if (!is.finite(curratio) || curratio==0) {
                stop("median ratio of averages between batches is not finite")
            }

            collected.ratios[first,second] <- curratio
            collected.ratios[second,first] <- 1/curratio
        }
    }

    # Finding the smallest ratio, and using that as the reference.
    smallest <- which.min(apply(collected.ratios, 2, min, na.rm=TRUE))
    
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
