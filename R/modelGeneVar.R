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

#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm 
.model_gene_var_per_block <- function(x, ..., design=NULL, subset.fit=NULL, subset.row=NULL, BPPARAM=SerialParam()) {
    fit <- fitTrendVar(x, ..., design=design, subset.row=subset.fit, BPPARAM=BPPARAM)

    if (identical(subset.fit, subset.row)) {
        x.mean <- fit$mean
        x.vars <- fit$var
    } else {
        out <- .get_var_stats(x, block=NULL, design=design, subset.row=subset.row, BPPARAM=BPPARAM)
        x.mean <- out$mean
        x.vars <- out$var
    }

    output <- DataFrame(mean=x.mean, total=x.vars, tech=fit$trend(x.mean))
    output$bio <- output$total - output$tech
    output$p.value <- pnorm(output$bio/output$tech, sd=fit$rsd, lower.tail=FALSE)

    rownames(output) <- rownames(x)[.subset_to_index(subset.row, x)]
    metadata(output) <- fit

    output
}

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pnorm p.adjust
.model_gene_var <- function(x, ..., design=NULL, subset.fit=NULL, subset.row=NULL, 
    block=NULL, equiweight=TRUE, method="fisher", BPPARAM=SerialParam()) 
{
    args <- list(..., subset.fit=subset.fit, subset.row=subset.row, BPPARAM=BPPARAM)

    if (is.null(block)) {
        output <- do.call(.model_gene_var_per_block, c(list(x=x, design=design), args))
        output$FDR <- p.adjust(output$p.value, method="BH")
        return(output)
    }

    collected <- split(seq_along(block), block)
    if (equiweight) {
        weights <- rep(1, length(collected))
    } else {
        weights <- lengths(collected)
    }

    for (i in names(collected)) {
        current <- collected[[i]]
        cur.args <- list(x=x[,current], design=design[current,,drop=FALSE])
        collected[[i]] <- do.call(.model_gene_var_per_block, c(cur.args, args))
    }

    # Combining statistics with optional weighting.
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
#' @rdname modelGeneVar
setMethod("modelGeneVar", "SingleCellExperiment", function(x, ..., subset.fit=NULL, assay.type="logcounts", use.spikes=TRUE) {
    if (use.spikes) {
        subset.fit <- isSpike(x)
        if (!any(subset.fit)) {
            stop("no spike-ins to use with 'use.spikes=TRUE'")
        }
    }
    .check_centered_SF(x, assay.type=assay.type)
    .model_gene_var(assay(x, i=assay.type), ..., subset.fit=subset.fit)
}) 
