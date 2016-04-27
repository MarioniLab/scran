.normalize.SCData <- function(object, size.factor=NULL, log=TRUE, prior.count=1) 
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
}

setMethod("normalize", "SCESet", function(object, log=TRUE, prior.count=1, separate.spikes=TRUE) {
    out <- .normalize.SCData(assayDataElement(object, "counts"), size.factor=sizeFactors(object), 
                             log=log, prior.count=prior.count) # Normalizing everything, not just spikes.

    if (separate.spikes && any(is.spike(object))) {
        object2 <- computeSpikeFactors(object)
        out2 <- .normalize.SCData(spikes(object, assay="counts"), size.factor=sizeFactors(object2), 
                                  log=log, prior.count=prior.count)
        out[is.spike(object),] <- out2
    } 
    
    assayDataElement(object, "exprs") <- out
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
