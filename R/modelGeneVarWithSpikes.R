#' Model the per-gene variance with spike-ins
#'
#' Model the variance of the log-expression profiles for each gene, 
#' decomposing it into technical and biological components based on a mean-variance trend fitted to spike-in transcripts.
#' 
#' @param x A numeric matrix of counts where rows are (usually endogenous) genes and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param spikes A numeric matrix of counts where each row corresponds to a spike-in transcript.
#' This should have the same number of columns as \code{x}.
#'
#' Alternatively, for the SingleCellExperiment method, 
#' this can be a string or an integer scalar specifying the \code{\link{altExp}} containing the spike-in count matrix.
#' @param size.factors A numeric vector of size factors for each cell in \code{x}, to be used for scaling gene expression.
#' @param spike.size.factors A numeric vector of size factors for each cell in \code{spikes}, to be used for scaling spike-ins.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add prior to log-transformation.
#' @param ... For the generic, further arguments to pass to each method.
#' 
#' For the ANY method, further arguments to pass to \code{\link{fitTrendVar}}.
#'
#' For the \linkS4class{SummarizedExperiment} method, further arguments to pass to the ANY method.
#'
#' For the \linkS4class{SingleCellExperiment} method, further arguments to pass to the SummarizedExperiment method.
#' @param assay.type String or integer scalar specifying the assay containing the counts.
#'
#' For the SingleCellExperiment method, this is used to retrieve both the endogenous and spike-in counts.
#' @inheritParams modelGeneVar
#'
#' @details
#' For each gene and spike-in transcript, we compute the variance and mean of the log-expression values.
#' A trend is fitted to the variance against the mean for spike-in transcripts using \code{\link{fitTrendVar}}.
#' The technical component for each gene is defined as the value of the trend at that gene's mean abundance.
#' The biological component is then defined as the residual from the trend.
#'
#' This function is almost the same as \code{\link{modelGeneVar}}, with the only difference being that the trend is fitted on spike-in variances rather than using the means and variances of endogenous genes.
#' It assumes that a constant amount of spike-in RNA was added to each cell, such that any differences in observed expression of the spike-in transcripts can be wholly attributed to technical noise; 
#' and that endogenous genes and spike-in transcripts follow the same mean-variance relationship.
#' 
#' Unlike \code{\link{modelGeneVar}}, \code{modelGeneVarWithSpikes} starts from a count matrix (for both genes and spike-ins) and requires size factors and a pseudo-count specification to compute the log-expression values.
#' This is because there are certain requirements for how normalization is performed when comparing spike-in transcripts with endogenous genes - see comments in \dQuote{Explanation for size factor rescaling}.
#' We enforce this by centering the size factors for both sets of features and recomputing the log-expression values prior to computing means and variances.
#'
#' @section Default size factor choices:
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
#' @section Explanation for size factor rescaling:
#' The use of a spike-in-derived trend makes several assumptions.
#' The first is that a constant amount of spike-in RNA was added to each cell, such that any differences in observed expression of the spike-in transcripts can be wholly attributed to technical noise. 
#' The second is that endogenous genes and spike-in transcripts follow the same mean-variance relationship, i.e., a spike-in transcript captures the technical noise of an endogenous gene at the same mean count.
#' 
#' Here, the spike-in size factors across all cells are scaled so that their mean is equal to that of the gene-based size factors for the same set of cells.
#' This ensures that the average normalized abundances of the spike-in transcripts are comparable to those of the endogenous genes,
#' allowing the trend fitted to the former to be used to determine the biological component of the latter. 
#' Otherwise, differences in scaling of the size factors would shift the normalized expression values of the former away from the latter, violating the second assumption.
#'
#' If \code{block} is specified, rescaling is performed separately for all cells in each block.
#' This aims to avoid problems from (frequent) violations of the first assumption where there are differences in the quantity of spike-in RNA added to each batch.
#' Without scaling, these differences would lead to systematic shifts in the spike-in abundances from the endogenous genes when fitting a batch-specific trend
#' (even if there is no global difference in scaling across all batches).
#' 
#' @inheritSection modelGeneVar Handling uninteresting factors
#' 
#' @section Computing p-values:
#' The p-value for each gene is computed by assuming that the variance estimates are normally distributed around the trend, and that the standard deviation of the variance distribution is proportional to the value of the trend.
#' This is used to construct a one-sided test for each gene based on its \code{bio}, under the null hypothesis that the biological component is equal to zero.
#' The proportionality constant for the standard deviation is set to the \code{std.dev} returned by \code{\link{fitTrendVar}}.
#' This is estimated from the spread of variance estimates for spike-in transcripts, so the null hypothesis effectively becomes \dQuote{is this gene \emph{more} variable than spike-in transcripts of the same abundance?}
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
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendVar}} on the spike-in transcripts,
#' along with the \code{mean} and \code{var} used to fit the trend.
#'
#' If \code{block} is specified, the output contains another \code{per.block} field.
#' This field is itself a DataFrame of DataFrames, where each internal DataFrame contains statistics for the variance modelling within each block and has the same format as described above. 
#' Each internal DataFrame's \code{metadata} contains the output of \code{\link{fitTrendVar}} for the cells of that block.
#'
#' @author Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Using spike-ins.
#' spk <- modelGeneVarWithSpikes(sce, "Spikes")
#' spk
#' 
#' plot(spk$mean, spk$total, log="xy")
#' points(metadata(spk)$mean, metadata(spk)$cv2, col="red", pch=16)
#' curve(metadata(spk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking (and spike-ins).
#' block <- sample(LETTERS[1:2], ncol(sce), replace=TRUE)
#' blk <- modelGeneVarWithSpikes(sce, "Spikes", block=block)
#' blk
#'
#' par(mfrow=c(1,2))
#' for (i in colnames(blk$per.block)) {
#'     current <- blk$per.block[[i]]
#'     plot(current$mean, current$total)
#'     points(metadata(current)$mean, metadata(current)$cv2, col="red", pch=16)
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
#' @importFrom scuttle librarySizeFactors .subset2index
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

    collected <- .decompose_log_exprs(all$x$means, all$x$vars, all$spikes$means, all$spikes$vars, 
        all$x$ncells, ...)
    output <- .combine_blocked_statistics(collected, method, equiweight, all$x$ncells)
    rownames(output) <- rownames(x)[.subset2index(subset.row, x)]
    output
}

#########################
# Setting up S4 methods #
#########################

#' @export
#' @rdname modelGeneVarWithSpikes
setGeneric("modelGeneVarWithSpikes", function(x, ...) standardGeneric("modelGeneVarWithSpikes"))

#' @export
#' @rdname modelGeneVarWithSpikes
setMethod("modelGeneVarWithSpikes", "ANY", .model_gene_var_with_spikes)

#' @export
#' @importFrom SummarizedExperiment assay
#' @rdname modelGeneVarWithSpikes
setMethod("modelGeneVarWithSpikes", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .model_gene_var_with_spikes(x=assay(x, i=assay.type), ...)
})

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

    callNextMethod(x=x, assay.type=assay.type, spikes=spikes,
        size.factors=size.factors, spike.size.factors=spike.size.factors, ...)
}) 
