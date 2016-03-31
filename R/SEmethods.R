setMethod("normalize", "ANY", function(object, size.factor=NULL, log=TRUE, prior.count=1) 
# Computes the normalized log-expression values.
# 
# written by Aaron Lun
# created 17 February 2016
# last modified 15 March 2016
{
    object <- as.matrix(object)
    if (is.null(size.factor)) { size.factor <- colSums(object) } 
    lsf <- log(size.factor) # Mean-centered size factors, for valid comparisons between size factor sets.
    size.factor <- exp(lsf - mean(lsf))
    cpm.default(object, lib.size=size.factor*1e6, prior.count=prior.count, log=log)
})

setMethod("normalize", "SCESet", function(object, ..., separate.spikes=TRUE) {
    out <- normalize(assayDataElement(object, "counts"), size.factor=sizeFactors(object), ...) # Normalizing everything, not just spikes.

    if (separate.spikes && any(is.spike(object))) {
        object2 <- computeSpikeFactors(object)
        out2 <- normalize(spikes(object, assay="counts"), size.factor=sizeFactors(object2), ...)
        out[is.spike(object),] <- out2
    } 
    
    assayDataElement(object, "exprs") <- out
    return(object)
})

setMethod("sizeFactors", "SCESet", function(object) {
    out <- object$sizeFactor
    if (is.null(out)) { 
        warning("'sizeFactors' are not set, returning NULL") 
        return(NULL)
    }
    names(out) <- colnames(object) 
    return(out)
})

setReplaceMethod("sizeFactors", "SCESet", function(object, value) {
    if (!is.numeric(value) && !is.null(value)) { 
        stop("size factors should be numeric or NULL")
    }        
    object$sizeFactor <- value
    return(object)
})

.getUsedMatrix <- function(x, assay="counts", get.spikes=FALSE) {
    cur.mat <- assayDataElement(x, assay)
    if (!get.spikes) {
        nokeep <- is.spike(x)
        if (!is.null(nokeep) && any(nokeep)) { 
            cur.mat <- cur.mat[!nokeep,,drop=FALSE]
        }
    }
    return(cur.mat)
}
