#' @export
setGeneric("computeSpikeFactors", function(x, ...) { standardGeneric("computeSpikeFactors") })

#' @importFrom SingleCellExperiment isSpike spikeNames
#' @importFrom BiocGenerics sizeFactors "sizeFactors<-"
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2
#' @export
setMethod("computeSpikeFactors", "SingleCellExperiment", 
          function(x, type=NULL, assay.type="counts", sf.out=FALSE, general.use=TRUE) 
# Uses the mean-centred total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
{
    .Deprecated(msg="'computeSpikeFactors' is not necessary if spike-ins are stored as 'altExps'")
    if (is.null(type)) { 
        is.spike <- isSpike(x)
    } else {
        is.spike <- logical(nrow(x))
        for (tset in type) {
            current <- isSpike(x, type=tset)
            if (!is.null(current)) { is.spike <- is.spike | current }
        }
    }
    if (!any(is.spike)) {
        stop("no spike-in transcripts present in 'x'")
    }

    # Computing spike-in size factors.
    out <- colSums2(DelayedArray(assay(x, i=assay.type)), rows=is.spike)
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
        type <- spikeNames(x)
    }
    for (f in type) {
        sizeFactors(x, type=f) <- sf
    }        
    x
})

