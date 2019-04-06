#' @export
#' @importFrom BiocGenerics normalize
#' @importFrom SingleCellExperiment isSpike 
multiBatchNorm <- function(..., assay.type="counts", norm.args=list(), min.mean=1, subset.row=NULL)
# Performs multi-batch normalization, adjusting for differences 
# in scale between SCE objects supplied in '...'.
# 
# written by Aaron Lun
# created 4 June 2018
{
    .Deprecated("batchelor::multiBatchNorm")
    batches <- list(...)
    .check_batch_consistency(batches, byrow=TRUE)
    .check_spike_consistency(batches)
    if (length(batches)==0L) {
        stop("at least one SingleCellExperiment object must be supplied")
    }

    # Adjusting size factors for the non-spike-in transcripts.
    nonspike.subset <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes=FALSE)
    if (is.null(nonspike.subset) || any(nonspike.subset)) { # works logical or integer
        batches <- .batch_rescaler(batches, subset.row=nonspike.subset, exprs_values=assay.type, sf.type=NULL, min.mean=min.mean)
    }

    # Adjusting size factors for each of the spike-ins.
    ref.spike.names <- spikeNames(batches[[1]])
    subset.row <- .subset_to_index(subset.row, batches[[1]], byrow=TRUE)

    for (spike in ref.spike.names) {
        ref.spike <- isSpike(batches[[1]], spike)
        for (b in seq_along(batches)) {
            if (!identical(ref.spike, isSpike(batches[[b]], spike))) {
                stop(sprintf("%s spike-in identities differ across batches", spike))
            }
        }

        spike.subset <- intersect(which(ref.spike), subset.row)
        if (length(spike.subset)) { 
            batches <- .batch_rescaler(batches, subset.row=spike.subset, exprs_values=assay.type, sf.type=spike, min.mean=min.mean)
        }
    }

    # Applying the normalization.   
    mapply(FUN=normalize, batches, MoreArgs=c(list(exprs_values=assay.type, centre_size_factors=FALSE), norm.args), SIMPLIFY=FALSE)
}

#' @importFrom scater calcAverage
#' @importFrom stats median
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
.batch_rescaler <- function(batches, subset.row, exprs_values, sf.type, min.mean) 
# Computes the median ratios (a la DESeq normalization),
# finds the smallest ratio and uses that as the reference.
{
    collected.ave <- lapply(batches, calcAverage, use_size_factors=FALSE, subset_row=subset.row, exprs_values=exprs_values)
    nbatches <- length(batches)
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

    # Rescaling to the lowest coverage.
    smallest <- which.min(apply(collected.ratios, 2, min, na.rm=TRUE))
    
    for (idx in seq_along(batches)) {
        current <- batches[[idx]]
        cursf <- sizeFactors(current, sf.type)
        if (is.null(cursf)) {
            cursf <- librarySizeFactors(current, exprs_values=exprs_values, subset_row=subset.row)
            warning(sprintf("no %ssize factors in batch %i, using sum of counts instead",
                if (is.null(sf.type)) "" else paste(sf.type, ""), idx))
        }
    
        # Adjusting the size factors for the new reference.
        cursf <- cursf/mean(cursf) / collected.ratios[idx,smallest]
        sizeFactors(current, sf.type) <- cursf
        batches[[idx]] <- current
    }

    return(batches)
}
