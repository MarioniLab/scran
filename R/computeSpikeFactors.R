#' Normalization with spike-in counts
#'
#' Compute size factors based on the coverage of spike-in transcripts.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object containing spike-in transcripts in an \code{\link{altExps}} entry.
#'
#' Support for spike-in transcripts in the rows of \code{x} itself is deprecated.
#' @param spikes String or integer scalar specifying the alternative experiment containing the spike-in transcripts.
#' @param type Deprecated, a character vector specifying which spike-in sets to use.
#' @param assay.type A string indicating which assay contains the counts.
#' @param sf.out Deprecated, a logical scalar indicating whether only size factors should be returned.
#' @param general.use A logical scalar indicating whether the size factors should be stored for general use by all genes.
#' 
#' @details
#' The spike-in size factor for each cell is computed from the sum of all spike-in counts in each cell.
#' This aims to scale the counts to equalize spike-in coverage between cells,
#' thus removing differences in coverage due to technical effects like capture or amplification efficiency.
#' 
#' Spike-in normalization can be helpful for preserving changes in total RNA content between cells, if this is of interest.
#' Such changes would otherwise be lost when normalizing with methods that assume a non-DE majority.
#' Indeed, spike-in normalization is the only available approach if a majority of genes are DE between two cell types or states.
#'
#' Size factors are computed by applying \code{\link{librarySizeFactors}} to the spike-in count matrix.
#' This ensures that the mean of all size factors is unity for standardization purposes,
#' if one were to compare expression values normalized with sets of size factors (e.g., in \code{\link{modelGeneVarWithSpikes}}).
#'
#' @return
#' A modified \code{x} is returned, containing spike-in-derived size factors for all cells in \code{\link{sizeFactors}}.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- computeSpikeFactors(sce, "Spikes")
#' summary(sizeFactors(sce))
#' 
#' # Deprecated behavior:
#' y <- sce
#' altExps(y) <- NULL
#' y <- rbind(y, altExp(sce))
#' suppressWarnings(isSpike(y, "ERCC") <- nrow(sce) + seq_len(nrow(altExp(sce))))
#' suppressWarnings(y <- computeSpikeFactors(y))
#' 
#' y2 <- y
#' @seealso
#' \code{\link{altExps}}, for the concept of alternative experiments.
#'
#' \code{\link{librarySizeFactors}}, for how size factors are derived from library sizes.
#' 
#' @references
#' Lun ATL, McCarthy DJ and Marioni JC (2016). 
#' A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.
#' \emph{F1000Res.} 5:2122
#' 
#' @export
#' @importFrom SingleCellExperiment isSpike spikeNames altExp
#' @importFrom BiocGenerics sizeFactors "sizeFactors<-"
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2
computeSpikeFactors <- function(x, spikes, type=NULL, assay.type="counts", sf.out=FALSE, general.use=TRUE) {
    if (!missing(spikes)) {
        sizeFactors(x) <- librarySizeFactors(altExp(x, spikes), exprs_values=assay.type)
        return(x)
    }

    if (is.null(type)) { 
        suppressWarnings(is.spike <- isSpike(x))
    } else {
        is.spike <- logical(nrow(x))
        for (tset in type) {
            suppressWarnings(current <- isSpike(x, type=tset))
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
        .Deprecated(msg="'general.use=TRUE' is deprecated.\nUse 'spikes=' instead, for spike-ins stored as 'altExps'.")
    } else {
        .Deprecated(msg="'computeSpikeFactors' is not necessary if spike-ins are stored as 'altExps'")
    }
    if (is.null(type)) {
        suppressWarnings(type <- spikeNames(x))
    }
    for (f in type) {
        suppressWarnings(sizeFactors(x, type=f) <- sf)
    }        
    x
}
