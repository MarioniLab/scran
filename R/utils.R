.get_feature_control_names <- function(x) {
    x@featureControlInfo$name
}

.get_feature_control_spike_names <- function(x) {
    x@featureControlInfo$name[x@featureControlInfo$spike]
}

is.spike <- function(x, type=NULL) { 
    if (is.null(type)) {
        out <- fData(x)$is_feature_spike
        return(out)
    } else {
        not.in <- !(type %in% .get_feature_control_spike_names(x))
        if (any(not.in)) { 
            stop(sprintf("'%s' is not specified as a spike-in control", type[which(not.in)[1]]))
        } 
        
        # Returning directly if possible.
        if (length(type)==1L) {
            return(fData(x)[[paste0("is_feature_control_", type)]])
        }

        # Combining the spike-in identities. 
        is.spike <- logical(nrow(x)) 
        for (f in type) {
            is.spike <- is.spike | fData(x)[[paste0("is_feature_control_", f)]]
        }
        return(is.spike)
    }
}

.spikeSubset <- function(x, get.spikes) {
    if (!get.spikes) {
        nokeep <- is.spike(x)
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
