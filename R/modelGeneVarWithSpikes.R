#' Model the per-gene variance
#'
#' Model the variance of the log-expression profiles for each gene, 
#' decomposing it into technical and biological components based on a fitted mean-variance trend.
#' 
#' @param x A numeric matrix of log-counts, or a \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param design A numeric matrix containing blocking terms for uninteresting factors of variation.
#' @param min.mean A numeric scalar specifying the minimum mean to use for trend fitting.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}, specifying the rows for which to model the variance.
#' @param subset.fit An argument similar to \code{subset.row}, specifying the rows to be used for trend fitting.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be performed across genes.
#' @param ... For the \code{ANY} method, further arguments to pass to \code{\link{fitTrendVar}}.
#'
#' For the generic and \linkS4class{SingleCellExperiment} method, further arguments to pass to the \code{ANY} method.
#' @param block A factor specifying the blocking levels for each cell in \code{x}.
#' If specified, variance modelling is performed separately within each block and statistics are combined across blocks.
#' @param equiweight A logical scalar indicating whether statistics from each block should be given equal weight.
#' Otherwise, each block is weighted according to its number of cells.
#' Only used if \code{block} is specified.
#' @param method String specifying how p-values should be combined, see \code{\link{combinePValues}}.
#' @param assay.type String or integer scalar specifying the assay containing the log-expression values.
#' @param use.spikes Logical scalar indicating whether spike-in transcripts should be used for trend fitting.
#'
#' @details
#' For each gene, the technical component of the variance is estimated according to the fitted value of the trend.
#' The biological component is defined as the residual from the trend.
#' This enables the identification of interesting genes in a manner that accounts for the mean-variance relationship.
#'
#' @return 
#' A \linkS4class{DataFrame} is returned where each row corresponds to a gene in \code{x} (or in \code{subset.row}, if specified).
#' This contains the numeric fields:
#' \describe{
#' \item{\code{mean}:}{Mean normalized log-expression per gene.}
#' \item{\code{total}:}{Variance of the normalized log-expression per gene.}
#' \item{\code{bio}:}{Biological component of the variance.}
#' \item{\code{tech}:}{Technical component of the variance.}
#' \item{\code{p.value, FDR}:}{Raw and adjusted p-values for the test against the null hypothesis that \code{bio=0}.}
#' }
#' 
#' If \code{block} is not specified, 
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendVar}} on the specified features.
#'
#' If \code{block} is specified,
#' the output contains another \code{per.block} field.
#' This field is itself a DataFrame of DataFrames, where each internal DataFrame contains statistics for the variance modelling within each block and has the same format as described above. 
#' Each internal DataFrame's \code{metadata} contains the output of \code{\link{fitTrendVar}} for the cells of that block.
#'
#' @author Aaron Lun
#' 
#' @examples
#' data(example.sce)
#'
#' # Using spike-ins.
#' spk <- modelGeneVar(example.sce)
#' spk
#' 
#' plot(spk$mean, spk$total)
#' points(metadata(spk)$mean, metadata(spk)$var, col="red")
#' curve(metadata(spk)$trend(x), add=TRUE, col="blue")
#'
#' # Not using spike-ins.
#' nspk <- modelGeneVar(example.sce, use.spikes=FALSE)
#' nspk
#' 
#' plot(nspk$mean, nspk$total)
#' curve(metadata(nspk)$trend(x), add=TRUE, col="blue")
#'
#' # With blocking (and spike-ins).
#' block <- sample(LETTERS[1:2], ncol(example.sce), replace=TRUE)
#' blk <- modelGeneVar(example.sce, block=block)
#' blk
#'
#' par(mfrow=c(1,2))
#' for (i in colnames(blk$per.block)) {
#'     current <- blk$per.block[[i]]
#'     plot(current$mean, current$total)
#'     points(metadata(current)$mean, metadata(current)$var, col="red")
#'     curve(metadata(current)$trend(x), add=TRUE, col="blue")
#' }
#' 
#' @name modelGeneVarWithSpikes
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
