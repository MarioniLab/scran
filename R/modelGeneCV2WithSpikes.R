#' Model the per-gene CV2 with spike-ins
#'
#' Model the squared coefficient of variation (CV2) of the normalized expression profiles for each gene, 
#' using spike-ins to estimate the baseline level of technical noise at each abundance.
#' 
#' @param ... For the generic, further arguments to pass to each method.
#' 
#' For the ANY method, further arguments to pass to \code{\link{fitTrendCV2}}.
#'
#' For the \linkS4class{SummarizedExperiment} method, further arguments to pass to the ANY method.
#'
#' For the \linkS4class{SingleCellExperiment} method, further arguments to pass to the SummarizedExperiment method.
#' @inheritParams modelGeneVarWithSpikes
#'
#' @details
#' For each gene and spike-in transcript, we compute the variance and CV2 of the normalized expression values.
#' A trend is fitted to the CV2 against the mean for spike-in transcripts using \code{\link{fitTrendCV2}}.
#' The value of the trend at the abundance of each gene is used to define the variation attributable to technical noise.
#' The ratio to the trend is then used to define overdispersion corresponding to interesting biological heterogeneity.
#' 
#' This function is almost the same as \code{\link{modelGeneCV2}}, with the only theoretical difference being that the trend is fitted on spike-in CV2 rather than using the means and CV2 of endogenous genes.
#' This is because there are certain requirements for how normalization is performed when comparing spike-in transcripts with endogenous genes - see comments in \dQuote{Explanation for size factor rescaling}.
#' We enforce this by centering the size factors for both sets of features and recomputing normalized expression values. 
#'
#' @inheritSection modelGeneVarWithSpikes Default size factor choices
#'
#' @inheritSection modelGeneVarWithSpikes Explanation for size factor rescaling
#'
#' @inheritSection modelGeneCV2 Handling uninteresting factors
#'
#' @section Computing p-values:
#' The p-value for each gene is computed by assuming that the CV2 estimates are normally distributed around the trend, and that the standard deviation of the CV2 distribution is proportional to the value of the trend.
#' This is used to construct a one-sided test for each gene based on its \code{ratio}, under the null hypothesis that the ratio is equal to 1.
#' The proportionality constant for the standard deviation is set to the \code{std.dev} returned by \code{\link{fitTrendCV2}}.
#' This is estimated from the spread of CV2 values for spike-in transcripts, so the null hypothesis effectively becomes \dQuote{is this gene \emph{more} variable than spike-in transcripts of the same abundance?}
#' 
#' @return 
#' A \linkS4class{DataFrame} is returned where each row corresponds to a gene in \code{x} (or in \code{subset.row}, if specified).
#' This contains the numeric fields:
#' \describe{
#' \item{\code{mean}:}{Mean normalized expression per gene.}
#' \item{\code{total}:}{Squared coefficient of variation of the normalized expression per gene.}
#' \item{\code{trend}:}{Fitted value of the trend.}
#' \item{\code{ratio}:}{Ratio of \code{total} to \code{trend}.}
#' \item{\code{p.value, FDR}:}{Raw and adjusted p-values for the test against the null hypothesis that \code{ratio<=1}.}
#' }
#' 
#' If \code{block} is not specified, 
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendCV2}} on the spike-in transcripts,
#' along with the \code{mean} and \code{cv2} used to fit the trend.
#'
#' If \code{block} is specified, the output contains another \code{per.block} field.
#' This field is itself a DataFrame of DataFrames, where each internal DataFrame contains statistics for the variance modelling within each block and has the same format as described above. 
#' Each internal DataFrame's \code{metadata} contains the output of \code{\link{fitTrendCV2}} for the cells of that block.
#'
#' @author Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#'
#' # Using spike-ins.
#' spk <- modelGeneCV2WithSpikes(sce, "Spikes")
#' spk
#' 
#' plot(spk$mean, spk$total)
#' points(metadata(spk)$mean, metadata(spk)$var, col="red", pch=16)
#' curve(metadata(spk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking (and spike-ins).
#' block <- sample(LETTERS[1:2], ncol(sce), replace=TRUE)
#' blk <- modelGeneCV2WithSpikes(sce, "Spikes", block=block)
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
#' @name modelGeneCV2WithSpikes
#' @seealso
#' \code{\link{fitTrendCV2}}, for the trend fitting options.
#'
#' \code{\link{modelGeneCV2}}, for modelling variance without spike-in controls.
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pnorm p.adjust
#' @importFrom scuttle librarySizeFactors .subset2index
.model_gene_cv2_with_spikes <- function(x, spikes, size.factors=NULL, spike.size.factors=NULL, 
    block=NULL, subset.row=NULL, ..., 
    equiweight=TRUE, method="fisher", BPPARAM=SerialParam()) 
{
    all <- .compute_var_stats_with_spikes(x=x, size.factors=size.factors, 
        subset.row=subset.row, block=block, 
        spikes=spikes, spike.size.factors=spike.size.factors, 
        BPPARAM=BPPARAM,
        block.FUN=compute_blocked_stats_norm, 
        design=NULL)

    collected <- .decompose_cv2(all$x$means, all$x$vars, all$spikes$means, all$spikes$vars, 
        ncells=all$x$ncells, ...)
    output <- .combine_blocked_statistics(collected, method, equiweight, all$x$ncells,
        geometric=TRUE, fields=c("mean", "total", "trend", "ratio"))
    rownames(output) <- rownames(x)[.subset2index(subset.row, x)]
    output
}

#########################
# Setting up S4 methods #
#########################

#' @export
#' @rdname modelGeneCV2WithSpikes
setGeneric("modelGeneCV2WithSpikes", function(x, ...) standardGeneric("modelGeneCV2WithSpikes"))

#' @export
#' @rdname modelGeneCV2WithSpikes
setMethod("modelGeneCV2WithSpikes", "ANY", .model_gene_cv2_with_spikes)

#' @export
#' @importFrom SummarizedExperiment assay
#' @rdname modelGeneCV2WithSpikes
setMethod("modelGeneCV2WithSpikes", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .model_gene_cv2_with_spikes(x=assay(x, i=assay.type), ...)
}) 

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment altExp
#' @importFrom methods selectMethod
#' @rdname modelGeneCV2WithSpikes
setMethod("modelGeneCV2WithSpikes", "SingleCellExperiment", function(x, spikes, 
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
