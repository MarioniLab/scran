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
#' @name modelGeneVar
#' @aliases modelGeneVar modelGeneVar,ANY-method modelGeneVar,SingleCellExperiment-method
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pnorm p.adjust
#' @importFrom scater librarySizeFactors
.model_gene_var <- function(x, size.factors=NULL, ..., design=NULL, 
    subset.row=NULL, block=NULL, fit.x=NULL, fit.size.factors=NULL, pseudo.count=1,
    equiweight=TRUE, method="fisher", BPPARAM=SerialParam()) 
{
    all <- .compute_var_stats_from_counts(x=x, size.factors=size.factors, 
        subset.row=subset.row, block=block, fit.x=fit.x, fit.size.factors=fit.size.factors, 
        FUN=.lognormvar, pseudo.count=pseudo.count, design=design)

    collected <- vector("list", ncol(all$x$means))
    for (i in seq_along(collected)) {
        fit <- .fit_trend_var0(all$fit$means[,i], all$fit$vars[,i], ...)

        x.means <- all$x$means[,i]
        output <- DataFrame(mean=x.means, total=all$x$vars[,i], tech=fit$trend(x.means))
        output$bio <- output$total - output$tech
        output$p.value <- pnorm(output$bio/output$tech, sd=fit$std.dev, lower.tail=FALSE)

        rownames(output) <- rownames(x.means)
        metadata(output) <- fit
        collected[[i]] <- output
    }

    if (length(collected)==1L) {
        return(collected[[1]])
    }

    # Combining statistics with optional weighting.
    if (equiweight) {
        weights <- rep(1, length(collected))
    } else {
        weights <- lengths(collected)
    }

    combined <- list()
    for (i in c("mean", "total", "tech", "bio")) {
        extracted <- lapply(collected, "[[", i=i)
        extracted <- mapply("*", extracted, weights, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        combined[[i]] <- Reduce("+", extracted)/sum(weights)
    }

    extracted <- lapply(collected, "[[", i="p.value")
    combined$p.value <- do.call(combinePValues, c(extracted, list(method=method, weights=weights)))
    combined$FDR <- p.adjust(combined$p.value, method="BH")

    output <- DataFrame(combined)
    rownames(output) <- rownames(x)[.subset_to_index(subset.row, x)]
    output$per.block <- do.call(DataFrame, lapply(collected, I))

    output
}

#########################
# Setting up S4 methods #
#########################

#' @export
setGeneric("modelGeneVar", function(x, ...) standardGeneric("modelGeneVar"))

#' @export
#' @rdname modelGeneVar
setMethod("modelGeneVar", "ANY", .model_gene_var)

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment altExp
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname modelGeneVar
setMethod("modelGeneVar", "SingleCellExperiment", function(x, size.factors=NULL, fit.x=NULL, fit.size.factors=NULL,
    ..., assay.type="counts", altexp=NULL) 
{
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    if (!is.null(altexp)){
        fit.x <- altExp(fit.x, alt.exp)
        if (is.null(fit.size.factors) && is(fit.x, "SingleCellExperiment")) {
            fit.size.factors <- sizeFactors(fit.x)            
        }
        fit.x <- assay(fit.x, assay.type)
    }
         
    .model_gene_var(x=assay(x, i=assay.type), size.factors=size.factors, 
        fit.x=fit.x, fit.size.factors=fit.size.factors, ...)
}) 
