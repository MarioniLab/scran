setGeneric("computeSpikeFactors", function(x, ...) { standardGeneric("computeSpikeFactors") })

setMethod("computeSpikeFactors", "SCESet", function(x, type=NULL, sf.out=FALSE, general.use=TRUE) 
# Uses the mean-centred total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
# last modified 28 May 2016
{
    out <- colSums(spikes(x, type=type))
    if (any(out < 1e-8)) { 
        warning("zero spike-in counts during spike-in normalization")
    } 
    sf <- out/mean(out)

    # Returning size factors directly.
    if (sf.out) {
        return(sf)
    }

    # Saving size factors for general use, or for specific use by one (or all) of the spike-in sets.
    if (general.use) {
        sizeFactors(x) <- sf
    } 
    if (is.null(type)) {
        type <- .get_feature_control_spike_names(x)
    }
    for (f in type) {
        sizeFactors(x, type=f) <- sf
    }        
    x
})

setGeneric("spikes", function(x, ...) standardGeneric("spikes"))

setMethod("spikes", "SCESet", function(x, assay="counts", type=NULL) {
    cur.assay <- assayDataElement(x, assay)[isSpike(x, type=type),,drop=FALSE]
    return(cur.assay)
})

setGeneric("isSpike", function(x, ...) standardGeneric("isSpike"))

setMethod("isSpike", "SCESet", function(x, type=NULL) {
    keep <- is.spike(x, type=type)
    if (is.null(keep)) {
        if (!is.null(type)) {
            extra <- sprintf(" for type='%s'", type)
        } else {
            extra <- ""
        }
        warning(sprintf("'isSpike' is not set%s, returning NULL", extra)) 
    }
    return(keep)
})

setGeneric("isSpike<-", function(x, value) standardGeneric("isSpike<-"))

setReplaceMethod("isSpike", signature(x="SCESet", value="NULL"), function(x, value) {
    fData(x)$is_feature_spike <- NULL 
    x@featureControlInfo$spike <- NULL
    return(x) 
})

setReplaceMethod("isSpike", signature(x="SCESet", value="character"), function(x, value) {
    # Recording all those that were listed as spikes.
    x@featureControlInfo$spike <- .get_feature_control_names(x) %in% value

    # Running through and collecting them.
    fData(x)$is_feature_spike <- is.spike(x, value)

    # Checking that they don't overlap.
    if (length(value) > 1L) { 
        total.hits <- integer(nrow(x))
        for (v in value) {
            total.hits <- total.hits + is.spike(x, v)
        }
        if (any(total.hits > 1L)) { 
            warning("overlapping spike-in sets detected")
        }
    }

    return(x) 
})

