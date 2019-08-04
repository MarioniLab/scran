#' Model the per-gene CV2 
#'
#' Model the squared coefficient of variation of the normalized expression profiles for each gene.
#' 
#' @param x A numeric matrix of counts, or a \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param size.factors A numeric vector of size factors for each cell/column in \code{x}.
#' @param fit.size.factors A numeric vector of size factors to apply to counts for features that are used in trend fitting.
#' This may be different from \code{size.factors} if spike-ins are to be used.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}, specifying the rows for which to model the variance.
#' @param subset.fit An argument similar to \code{subset.row}, specifying the rows to be used for trend fitting.
#' @param ... For the \code{ANY} method, further arguments to pass to \code{\link{fitTrendCV2}}.
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
#' For each gene, the mean-dependent trend in the CV2 values is modelled using \code{\link{fitTrendCV2}}.
#' Genes are considered to be more interesting if their CV2 values are much higher than the trend.
#' This is quantified based on the ratio of each gene's CV2 from the fitted value of the trend at the same abundance.
#'
#' @return 
#' A \linkS4class{DataFrame} is returned where each row corresponds to a gene in \code{x} (or in \code{subset.row}, if specified).
#' This contains the numeric fields:
#' \describe{
#' \item{\code{mean}:}{Mean normalized expression per gene.}
#' \item{\code{cv2}:}{CV2 of the normalized expression per gene.}
#' \item{\code{trend}:}{Fitted value of the trend.}
#' \item{\code{ratio}:}{Ratio of \code{cv2} to \code{trend}.}
#' \item{\code{p.value, FDR}:}{Raw and adjusted p-values for the test against the null hypothesis that \code{ratio=1}.}
#' }
#' 
#' If \code{block} is not specified, 
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendCV2}} on the specified features.
#'
#' If \code{block} is specified,
#' the output contains another \code{per.block} field.
#' This field is itself a DataFrame of DataFrames, where each internal DataFrame contains statistics for the variance modelling within each block and has the same format as described above. 
#' Each internal DataFrame's \code{metadata} contains the output of \code{\link{fitTrendCV2}} for the cells of that block.
#'
#' @author Aaron Lun
#' 
#' @examples
#' data(example.sce)
#'
#' # Simple case:
#' spk <- modelGeneCV2(example.sce)
#' spk
#' 
#' plot(spk$mean, spk$total, pch=16, log="xy")
#' curve(metadata(spk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking: 
#' block <- sample(LETTERS[1:2], ncol(example.sce), replace=TRUE)
#' blk <- modelGeneCV2(example.sce, block=block)
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
#' @importFrom scater librarySizeFactors
.model_gene_cv2 <- function(x, size.factors=NULL, block=NULL, subset.row=NULL, subset.fit=NULL,
    ..., equiweight=TRUE, method="fisher", BPPARAM=SerialParam())
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, subset_row=subset.row)
    }

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
    rownames(output) <- rownames(x)[.subset_to_index(subset.row, x)]
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
#' @importFrom BiocGenerics sizeFactors
#' @rdname modelGeneCV2
setMethod("modelGeneCV2", "SingleCellExperiment", function(x, size.factors=NULL, ..., assay.type="counts") {
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    .model_gene_cv2(assay(x, i=assay.type), ..., size.factors=size.factors)
}) 
