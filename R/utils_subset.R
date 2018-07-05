#' @importFrom SingleCellExperiment isSpike
.spike_subset <- function(x, get.spikes) 
# Returns a logical vector specifying which rows we should retain,
# if depending on whether or not we want to keep the spike-ins.
{
    if (!get.spikes) {
        nokeep <- isSpike(x)
        if (!is.null(nokeep) && any(nokeep)) {
            return(!nokeep)
        }
    } 
    return(NULL)
}

.subset_to_index <- function(subset, x, byrow=TRUE) 
# Converts arbitrary subsetting vectors to an integer index vector.
{
    if (byrow) {
        dummy <- seq_len(nrow(x))
        names(dummy) <- rownames(x)
    } else {
        dummy <- seq_len(ncol(x))
        names(dummy) <- colnames(x) 
    }

    if (!is.null(subset)) { 
        dummy <- dummy[subset]
    }
    out <- unname(dummy)
    if (any(is.na(out))) {
        stop("'subset' indices out of range of 'x'")
    }
    return(out)
}

.SCE_subset_genes <- function(subset.row, x, get.spikes) 
# Convenience function to intersect arbitrary subsetting specification with spike-in information.
# To be mainly used for SingleCellExperiments to further restrict subset.row.
{
    despiked <- .spike_subset(x, get.spikes)
    if (is.null(subset.row)) { 
        subset.row <- despiked
    } else {
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        if (!is.null(despiked)) { 
            subset.row <- intersect(subset.row, which(despiked))
        }
    }
    return(subset.row)
}
