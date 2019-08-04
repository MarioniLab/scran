#' Model the per-gene variance with spike-ins
#'
#' Model the variance of the log-expression profiles for each gene, 
#' decomposing it into technical and biological components based on a mean-variance trend fitted to spike-in transcripts.
#' 
#' @param x A numeric matrix of counts for endogenous genes, or a \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param spikes A numeric matrix of counts for spike-in transcripts.
#'
#' Alternatively, for the SingleCellExperiment method, this can be a string or an integer scalar specifying the \code{\link{altExp}} containing the spike-in count matrix.
#' @param size.factors A numeric vector of size factors for each cell in \code{x}, to be used for scaling gene expression.
#' @param spike.size.factors A numeric vector of size factors for each cell in \code{spikes}, to be used for scaling spike-ins.
#' @param design A numeric matrix containing blocking terms for uninteresting factors of variation.
#' @param min.mean A numeric scalar specifying the minimum mean to use for trend fitting.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}, specifying the rows for which to model the variance.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be performed across genes.
#' @param ... For the generic, further arguments to pass to each method.
#' 
#' For the ANY method, further arguments to pass to \code{\link{fitTrendVar}}.
#'
#' For the generic and \linkS4class{SingleCellExperiment} method, further arguments to pass to the ANY method.
#' @param block A factor specifying the blocking levels for each cell in \code{x}.
#' If specified, variance modelling is performed separately within each block and statistics are combined across blocks.
#' @param equiweight A logical scalar indicating whether statistics from each block should be given equal weight.
#' Otherwise, each block is weighted according to its number of cells.
#' Only used if \code{block} is specified.
#' @param method String specifying how p-values should be combined, see \code{\link{combinePValues}}.
#' @param assay.type String or integer scalar specifying the assay containing the counts.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add prior to log-transformation.
#'
#' @details
#' For each gene and spike-in transcript, we compute the variance and mean of the log-expression values.
#' A trend is fitted to the variance against the mean for spike-in transcripts using \code{\link{fitTrendVar}}.
#' The technical component for each gene is defined as the value of the trend at that gene's mean abundance,
#' under the assumption that both endogenous genes and spike-in transcripts follow the same mean-variance relationship.
#' The biological component is then defined as the residual from the trend.
#'
#' This function can be considered the same as \code{\link{modelGeneVar}}, with the only theoretical difference being that the trend is fitted on spike-in variances rather than using the means and variances of endogenous genes.
#' Both methods compute variances and means on the log-expression values;
#' both methods use \code{\link{fitTrendVar}} to fit the mean-variance trend; 
#' and both methods allow blocking or regression of uninteresting factors of variation.
#' 
#' Practically, \code{modelGeneVarWithSpikes} starts from a count matrix (for both genes and spike-ins) and requires size factors and a pseudo-count specification to compute the log-expression values.
#' In contrast, \code{\link{modelGeneVar}} simply starts from a log-expression matrix.
#' This discrepancy is deliberate and necessary to ensure that the log-expression values are comparable between genes and spike-ins.
#' Specifically, the mean size factor for the genes must be the same as the mean size factor for the spike-ins for the trend fitted to the latter to be applicable to the former.
#' (If \code{block} is specified, this must be true for all cells in each block.)
#' We enforce this by centering the size factors for both sets of features and recomputing the log-expression values.
#'
#' If no size factors are supplied, they are automatically computed depending on the input type:
#' \itemize{
#' \item If \code{size.factors=NULL} for the ANY method, the sum of counts for each cell in \code{x} is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item If \code{spike.size.factors=NULL} for the ANY method, the sum of counts for each cell in \code{spikes} is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item If \code{size.factors=NULL} for the \linkS4class{SingleCellExperiment} method, \code{\link{sizeFactors}(x)} is used if available.
#' Otherwise, it defaults to library size-derived size factors.
#' \item If \code{spike.size.factors=NULL} for the \linkS4class{SingleCellExperiment} method and \code{spikes} is not a matrix, \code{\link{sizeFactors}(\link{altExp}(x, spikes)} is used if available.
#' Otherwise, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} or \code{spike.size.factors} are supplied, they will override any size factors present in \code{x}.
#' 
#' @return 
#' A \linkS4class{DataFrame} is returned where each row corresponds to a gene in \code{x} (or in \code{subset.row}, if specified).
#' This contains the numeric fields:
#' \describe{
#' \item{\code{mean}:}{Mean normalized log-expression per gene.}
#' \item{\code{total}:}{Variance of the normalized log-expression per gene.}
#' \item{\code{bio}:}{Biological component of the variance.}
#' \item{\code{tech}:}{Technical component of the variance.}
#' \item{\code{p.value, FDR}:}{Raw and adjusted p-values for the test against the null hypothesis that \code{bio<=0}.}
#' }
#' 
#' If \code{block} is not specified, 
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendVar}} on the spike-in transcripts.
#'
#' If \code{block} is specified, the output contains another \code{per.block} field.
#' This field is itself a DataFrame of DataFrames, where each internal DataFrame contains statistics for the variance modelling within each block and has the same format as described above. 
#' Each internal DataFrame's \code{metadata} contains the output of \code{\link{fitTrendVar}} for the cells of that block.
#'
#' @author Aaron Lun
#' 
#' @examples
#' data(example.sce)
#'
#' # Using spike-ins.
#' spk <- modelGeneVarWithSpikes(example.sce, "Spike")
#' spk
#' 
#' plot(spk$mean, spk$total)
#' points(metadata(spk)$mean, metadata(spk)$var, col="red", pch=16)
#' curve(metadata(spk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking (and spike-ins).
#' block <- sample(LETTERS[1:2], ncol(example.sce), replace=TRUE)
#' blk <- modelGeneVarWithSpikes(example.sce, "Spike", block=block)
#' blk
#'
#' par(mfrow=c(1,2))
#' for (i in colnames(blk$per.block)) {
#'     current <- blk$per.block[[i]]
#'     plot(current$mean, current$total)
#'     points(metadata(current)$mean, metadata(current)$var, col="red", pch=16)
#'     curve(metadata(current)$trend(x), add=TRUE, col="dodgerblue")
#' }
#' 
#' @name modelGeneVarWithSpikes
#' @seealso
#' \code{\link{fitTrendVar}}, for the trend fitting options.
#'
#' \code{\link{modelGeneVar}}, for modelling variance without spike-in controls.
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pnorm p.adjust
#' @importFrom scater librarySizeFactors
.model_gene_var_with_spikes <- function(x, spikes, size.factors=NULL, spike.size.factors=NULL, 
    block=NULL, design=NULL, subset.row=NULL, pseudo.count=1, ..., 
    equiweight=TRUE, method="fisher", BPPARAM=SerialParam()) 
{
    all <- .compute_var_stats_with_spikes(x=x, size.factors=size.factors, 
        subset.row=subset.row, block=block, 
        spikes=spikes, spike.size.factors=spike.size.factors, 
        BPPARAM=BPPARAM,
        block.FUN=compute_blocked_stats_lognorm, 
        residual.FUN=compute_residual_stats_lognorm, 
        pseudo=pseudo.count, design=design)

    collected <- .decompose_log_exprs(all$x$means, all$x$vars, all$spikes$means, all$spikes$vars, ...)
    output <- .combine_blocked_statistics(collected, method, equiweight, all$x$ncells)
    rownames(output) <- rownames(x)[.subset_to_index(subset.row, x)]
    output
}

#########################
# Setting up S4 methods #
#########################

#' @export
setGeneric("modelGeneVarWithSpikes", function(x, ...) standardGeneric("modelGeneVarWithSpikes"))

#' @export
#' @rdname modelGeneVar
setMethod("modelGeneVarWithSpikes", "ANY", .model_gene_var_with_spikes)

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment altExp
#' @importFrom methods selectMethod
#' @rdname modelGeneVarWithSpikes
setMethod("modelGeneVarWithSpikes", "SingleCellExperiment", function(x, spikes, 
    size.factors=NULL, spike.size.factors=NULL, ..., assay.type="counts")
{
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }

    if (is.null(dim(spikes))) {
        spikes <- altExp(x, spikes)
        if (is.null(spike.size.factors)) {
            sfFUN <- selectMethod("sizeFactors", class(spikes), optional=TRUE)
            if (!is.null(sfFUN)) {
                spike.size.factors <- sfFUN(spikes)
            }
        }
        spikes <- assay(spikes, assay.type)
    }

    .model_gene_var_with_spikes(x=assay(x, i=assay.type), spikes=spikes,
        size.factors=size.factors, spike.size.factors=spike.size.factors, ...)
}) 
