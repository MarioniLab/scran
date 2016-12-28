.spikeSubset <- function(x, get.spikes) {
    if (!get.spikes) {
        nokeep <- isSpike(x, warning=FALSE)
        if (!is.null(nokeep) && any(nokeep)) {
            return(!nokeep)
        }
    } 
    return(NULL)
}

.subset_to_index <- function(subset, x, byrow=TRUE) {
    if (byrow) {
        dimlen <- nrow(x)
        names <- rownames(x)
    } else {
        dimlen <- ncol(x)
        names <- colnames(x)
    }

    if (is.logical(subset)) { 
        if (length(subset)!=dimlen) {
            stop("subset vector is longer than matrix dimensions") 
        }
        subset <- which(subset)
    } else if (is.character(subset)) {
        subset <- match(subset, names)
        if (any(is.na(subset))) { 
            stop("missing names in subset vector")
        }
    } else if (is.null(subset)) {
        subset <- seq_len(dimlen)
    } else if (is.numeric(subset)) {
        subset <- as.integer(subset)
        if (min(subset) < 1L || max(subset) > dimlen) {
            stop("subset indices out of range")
        }
    } else {
        stop("unrecognized type of subset vector")
    }

    return(subset)
}

.check_centered_SF <- function(x, assay) {
    # Checks if 'exprs' was requested, and if it could have been computed from counts,
    # If so, then it checks if the size factors are centered.
    if (assay=="exprs" && 
        !is.null(get_exprs(x, "counts", warning=FALSE)) && 
        !areSizeFactorsCentred(x)) {
        warning("size factors not centred, run 'normalize()' first")
    }
    return(NULL)
}

