#' Normalization with spike-in counts
#'
#' Compute size factors based on the coverage of spike-in transcripts.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object containing spike-in transcripts in an \code{\link{altExps}} entry.
#'
#' Support for spike-in transcripts in the rows of \code{x} itself is deprecated.
#' @param spikes String or integer scalar specifying the alternative experiment containing the spike-in transcripts.
#' @param assay.type A string indicating which assay contains the counts.
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
#' Users who want the spike-in size factors without returning a SingleCellExperiment object can simply call
#' \code{\link{librarySizeFactors}(\link{altExp}(x, spikes))}, which gives the same result.
#' 
#' @return
#' A modified \code{x} is returned, containing spike-in-derived size factors for all cells in \code{\link{sizeFactors}}.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- computeSpikeFactors(sce, "Spikes")
#' summary(sizeFactors(sce))
#' 
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
#' @importFrom SingleCellExperiment altExp
#' @importFrom BiocGenerics sizeFactors<-
#' @importFrom scuttle librarySizeFactors
computeSpikeFactors <- function(x, spikes, assay.type="counts") {
    sizeFactors(x) <- librarySizeFactors(altExp(x, spikes), exprs_values=assay.type)
    x
}
