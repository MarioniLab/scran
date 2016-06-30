setGeneric("computeSpikeFactors", function(x, ...) { standardGeneric("computeSpikeFactors") })

setMethod("computeSpikeFactors", "SCESet", function(x, type=NULL, sf.out=FALSE) 
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
    if (sf.out) {
        return(sf)
    }
    sizeFactors(x) <- sf
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

setReplaceMethod("isSpike", signature(x="SCESet", value="logical"), function(x, value) {
    fData(x)$is_feature_spike <- value
    return(x)
})

setReplaceMethod("isSpike", signature(x="SCESet", value="NULL"), function(x, value) {
    fData(x)$is_feature_spike <- NULL 
    return(x) 
})

setReplaceMethod("isSpike", signature(x="SCESet", value="character"), function(x, value) {
    # Checking that the controls exist.
    check.spikes <- ! (value %in% .get_feature_control_names(x))
    if (any(check.spikes)) { 
        stop(sprintf("need to specify '%s' as a control in calculateQCMetrics", value[which(check.spikes)[1]]))
    }

    # Running through and collecting them.
    is.spike <- logical(nrow(x)) 
    for (v in value) {
        is.spike <- is.spike | fData(x)[[paste0("is_feature_control_", v)]]
    }
    fData(x)$is_feature_spike <- is.spike
    return(x) 
})

