#' Model the per-gene variance with Poisson noise
#'
#' Model the variance of the log-expression profiles for each gene, 
#' decomposing it into technical and biological components based on a mean-variance trend corresponding to Poisson noise.
#' 
#' @param x A numeric matrix of counts where rows are (usually endogenous) genes and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param size.factors A numeric vector of size factors for each cell in \code{x}, to be used for scaling gene expression.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add prior to log-transformation.
#' @param ... For the generic, further arguments to pass to each method.
#' 
#' For the ANY method, further arguments to pass to \code{\link{fitTrendVar}}.
#'
#' For the \linkS4class{SummarizedExperiment} method, further arguments to pass to the ANY method.
#'
#' For the \linkS4class{SingleCellExperiment} method, further arguments to pass to the SummarizedExperiment method.
#' @inheritParams modelGeneVar
#' @inheritParams fitTrendPoisson
#' @param assay.type String or integer scalar specifying the assay containing the counts.
#'
#' @details
#' For each gene, we compute the variance and mean of the log-expression values.
#' A trend is fitted to the variance against the mean for simulated Poisson counts as described in \code{\link{fitTrendPoisson}}.
#' The technical component for each gene is defined as the value of the trend at that gene's mean abundance.
#' The biological component is then defined as the residual from the trend.
#'
#' This function is similar to \code{\link{modelGeneVarWithSpikes}}, with the only difference being that the trend is fitted on simulated Poisson count-derived variances rather than spike-ins.
#' The assumption is that the technical component is Poisson-distributed, or at least negative binomial-distributed with a known constant dispersion.
#' This is useful for UMI count data sets that do not have spike-ins and are too heterogeneous to assume that most genes exhibit negligible biological variability.
#'
#' If no size factors are supplied, they are automatically computed depending on the input type:
#' \itemize{
#' \item If \code{size.factors=NULL} for the ANY method, the sum of counts for each cell in \code{x} is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item If \code{size.factors=NULL} for the \linkS4class{SingleCellExperiment} method, \code{\link{sizeFactors}(x)} is used if available.
#' Otherwise, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} are supplied, they will override any size factors present in \code{x}.
#'
#' @inheritSection modelGeneVar Handling uninteresting factors
#' 
#' @section Computing p-values:
#' The p-value for each gene is computed by assuming that the variance estimates are normally distributed around the trend, and that the standard deviation of the variance distribution is proportional to the value of the trend.
#' This is used to construct a one-sided test for each gene based on its \code{bio}, under the null hypothesis that the biological component is equal to zero.
#' The proportionality constant for the standard deviation is set to the \code{std.dev} returned by \code{\link{fitTrendVar}}.
#' This is estimated from the spread of variance estimates for the simulated Poisson-distributed counts, so the null hypothesis effectively becomes \dQuote{is this gene \emph{more} variable than a hypothetical gene with only Poisson noise?}
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
#' the \code{metadata} of the DataFrame contains the output of running \code{\link{fitTrendVar}} on the simulated counts,
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
#'
#' # Using spike-ins.
#' pois <- modelGeneVarByPoisson(sce)
#' pois
#' 
#' plot(pois$mean, pois$total, ylim=c(0, 10))
#' points(metadata(pois)$mean, metadata(pois)$var, col="red", pch=16)
#' curve(metadata(pois)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking.
#' block <- sample(LETTERS[1:2], ncol(sce), replace=TRUE)
#' blk <- modelGeneVarByPoisson(sce, block=block)
#' blk
#'
#' par(mfrow=c(1,2))
#' for (i in colnames(blk$per.block)) {
#'     current <- blk$per.block[[i]]
#'     plot(current$mean, current$total, ylim=c(0, 10))
#'     points(metadata(current)$mean, metadata(current)$var, col="red", pch=16)
#'     curve(metadata(current)$trend(x), add=TRUE, col="dodgerblue")
#' }
#' 
#' @name modelGeneVarByPoisson
#' @seealso
#' \code{\link{fitTrendVar}}, for the trend fitting options.
#'
#' \code{\link{modelGeneVar}}, for modelling variance without spike-in controls.
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle librarySizeFactors .subset2index
.model_gene_var_by_poisson <- function(x, size.factors=NULL, 
    block=NULL, design=NULL, subset.row=NULL, npts=1000, dispersion=0, pseudo.count=1, ..., 
    equiweight=TRUE, method="fisher", BPPARAM=SerialParam()) 
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, subset_row=subset.row)
    }
    x.stats <- .compute_mean_var(x, block=block, design=design, subset.row=subset.row, 
        block.FUN=compute_blocked_stats_lognorm, 
        residual.FUN=compute_residual_stats_lognorm, 
        BPPARAM=BPPARAM, sf=size.factors, pseudo=pseudo.count)

    xlim <- 2^range(x.stats$means[x.stats$means > 0]) - pseudo.count
    sim.out <- .generate_poisson_values(xlim, size.factors, block=block, design=design,
        npts=npts, dispersion=dispersion, pseudo.count=pseudo.count, BPPARAM=BPPARAM)

    collected <- .decompose_log_exprs(x.stats$means, x.stats$vars, sim.out$means, sim.out$vars, 
        x.stats$ncells, ...)
    output <- .combine_blocked_statistics(collected, method, equiweight, x.stats$ncells)
    rownames(output) <- rownames(x)[.subset2index(subset.row, x)]
    output
}

#########################
# Setting up S4 methods #
#########################

#' @export
#' @rdname modelGeneVarByPoisson
setGeneric("modelGeneVarByPoisson", function(x, ...) standardGeneric("modelGeneVarByPoisson"))

#' @export
#' @rdname modelGeneVarByPoisson
setMethod("modelGeneVarByPoisson", "ANY", .model_gene_var_by_poisson)

#' @export
#' @importFrom SummarizedExperiment assay
#' @rdname modelGeneVarByPoisson
setMethod("modelGeneVarByPoisson", "SummarizedExperiment", function(x, ..., assay.type="counts")
{
    .model_gene_var_by_poisson(assay(x, i=assay.type), ...)
}) 

#' @export
#' @importFrom BiocGenerics sizeFactors
#' @importFrom methods selectMethod
#' @rdname modelGeneVarByPoisson
setMethod("modelGeneVarByPoisson", "SingleCellExperiment", function(x, size.factors=sizeFactors(x), ...)
{
    callNextMethod(x=x, size.factors=size.factors, ...)
}) 
