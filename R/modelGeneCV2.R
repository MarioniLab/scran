#' Model the per-gene CV2 
#'
#' Model the squared coefficient of variation (CV2) of the normalized expression profiles for each gene,
#' fitting a trend to account for the mean-variance relationship across genes.
#' 
#' @param x A numeric matrix of counts where rows are genes and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param size.factors A numeric vector of size factors for each cell in \code{x}.
#' @param ... For the generic, further arguments to pass to each method.
#' 
#' For the ANY method, further arguments to pass to \code{\link{fitTrendCV2}}.
#'
#' For the \linkS4class{SummarizedExperiment} method, further arguments to pass to the ANY method.
#'
#' For the \linkS4class{SingleCellExperiment} method, further arguments to pass to the SummarizedExperiment method.
#' @param assay.type String or integer scalar specifying the assay containing the counts.
#' @inheritParams modelGeneVar
#'
#' @details
#' For each gene, we compute the CV2 and mean of the counts after scaling them by \code{size.factors}.
#' A trend is fitted to the CV2 against the mean for all genes using \code{\link{fitTrendCV2}}.
#' The fitted value for each gene is used as a proxy for the technical noise, assuming that most genes exhibit a low baseline level of variation that is not biologically interesting.
#' The ratio of the total CV2 to the trend is used as a metric to rank interesting genes, with larger ratios being indicative of strong biological heterogeneity.
#'
#' By default, the trend is fitted using all of the genes in \code{x}.
#' If \code{subset.fit} is specified, the trend is fitted using only the specified subset,
#' and the values of \code{trend} for all other genes are determined by extrapolation or interpolation.
#' This could be used to perform the fit based on genes that are known to have low variance, thus weakening the assumption above.
#' Note that this does not refer to spike-in transcripts, which should be handled via \code{\link{modelGeneCV2WithSpikes}}.
#'
#' If no size factors are supplied, they are automatically computed depending on the input type:
#' \itemize{
#' \item If \code{size.factors=NULL} for the ANY method, the sum of counts for each cell in \code{x} is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item If \code{size.factors=NULL} for the \linkS4class{SingleCellExperiment} method, \code{\link{sizeFactors}(x)} is used if available.
#' Otherwise, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} are supplied, they will override any size factors present in \code{x}.
#'
#' @section Computing p-values:
#' The p-value for each gene is computed by assuming that the CV2 estimates are normally distributed around the trend, and that the standard deviation of the CV2 distribution is proportional to the value of the trend.
#' This is used to construct a one-sided test for each gene based on its \code{ratio}, under the null hypothesis that the ratio is equal to or less than 1.
#' The proportionality constant for the standard deviation is set to the \code{std.dev} returned by \code{\link{fitTrendCV2}}.
#' This is estimated from the spread of per-gene CV2 values around the trend, so the null hypothesis effectively becomes \dQuote{is this gene \emph{more} variable than other genes of the same abundance?}
#' 
#' @section Handling uninteresting factors:
#' Setting \code{block} will estimate the mean and variance of each gene for cells in each level of \code{block} separately.
#' The trend is fitted separately for each level, and the variance decomposition is also performed separately.
#' Per-level statistics are then combined to obtain a single value per gene:
#' \itemize{
#' \item For means and CV2 values, this is done by taking the geometric mean across blocking levels.
#' If \code{equiweight=FALSE}, a weighted average is used where the value for each level is weighted by the number of cells.
#' By default, all levels are equally weighted when combining statistics.
#' \item Per-level p-values are combined using \code{\link{combineParallelPValues}} according to \code{method}.
#' By default, Fisher's method is used to identify genes that are highly variable in any batch.
#' Whether or not this is responsive to \code{equiweight} depends on the chosen method.
#' \item Blocks with fewer than 2 cells are completely ignored and do not contribute to the combined mean, variance component or p-value.
#' }
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
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendCV2}} on the specified features,
#' along with the \code{mean} and \code{cv2} used to fit the trend.
#'
#' If \code{block} is specified,
#' the output contains another \code{per.block} field.
#' This field is itself a DataFrame of DataFrames, where each internal DataFrame contains statistics for the variance modelling within each block and has the same format as described above. 
#' Each internal DataFrame's \code{metadata} contains the output of \code{\link{fitTrendCV2}} for the cells of that block.
#'
#' @author Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#'
#' # Simple case:
#' spk <- modelGeneCV2(sce)
#' spk
#' 
#' plot(spk$mean, spk$total, pch=16, log="xy")
#' curve(metadata(spk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking: 
#' block <- sample(LETTERS[1:2], ncol(sce), replace=TRUE)
#' blk <- modelGeneCV2(sce, block=block)
#' blk
#'
#' par(mfrow=c(1,2))
#' for (i in colnames(blk$per.block)) {
#'     current <- blk$per.block[[i]]
#'     plot(current$mean, current$total, pch=16, log="xy")
#'     curve(metadata(current)$trend(x), add=TRUE, col="dodgerblue")
#' }
#' 
#' @name modelGeneCV2
#' @aliases modelGeneCV2 modelGeneCV2,ANY-method modelGeneCV2,SingleCellExperiment-method
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle librarySizeFactors .subset2index
.model_gene_cv2 <- function(x, size.factors=NULL, block=NULL, subset.row=NULL, subset.fit=NULL,
    ..., equiweight=TRUE, method="fisher", BPPARAM=SerialParam())
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, subset_row=subset.row)
    }
    size.factors <- size.factors/mean(size.factors)

    FUN <- function(s) {
        .compute_mean_var(x, block=block, design=NULL, subset.row=s,
            block.FUN=compute_blocked_stats_norm,
            BPPARAM=BPPARAM, sf=size.factors)
    }
    x.stats <- FUN(subset.row)

    if (is.null(subset.fit)) {
        fit.stats <- x.stats
    } else {
        # Yes, we could do this more efficiently by rolling up 'subset.fit'
        # into 'subset.row' for a single '.compute_mean_var' call... but I CBF'd.
        fit.stats <- FUN(subset.fit)
    }

    collected <- .decompose_cv2(x.stats$means, x.stats$vars, fit.stats$means, fit.stats$vars, 
        ncells=x.stats$ncells, ...)
    output <- .combine_blocked_statistics(collected, method, equiweight, x.stats$ncells, 
        geometric=TRUE, fields=c("mean", "total", "trend", "ratio"))
    rownames(output) <- rownames(x)[.subset2index(subset.row, x)]
    output
}

#########################
# Setting up S4 methods #
#########################

#' @export
#' @rdname modelGeneCV2
setGeneric("modelGeneCV2", function(x, ...) standardGeneric("modelGeneCV2"))

#' @export
#' @rdname modelGeneCV2
setMethod("modelGeneCV2", "ANY", .model_gene_cv2)

#' @export
#' @importFrom SummarizedExperiment assay
#' @rdname modelGeneCV2
setMethod("modelGeneCV2", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .model_gene_cv2(assay(x, i=assay.type), ...)
})

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
#' @rdname modelGeneCV2
setMethod("modelGeneCV2", "SingleCellExperiment", function(x, size.factors=NULL, ...) {
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    callNextMethod(x=x, ..., size.factors=size.factors)
}) 
