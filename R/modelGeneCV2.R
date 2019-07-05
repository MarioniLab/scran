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
#' # Using spike-ins.
#' spk <- modelGeneCV2(example.sce)
#' spk
#' 
#' plot(spk$mean, spk$cv2, pch=16, log="xy")
#' points(metadata(spk)$mean, metadata(spk)$cv2, col="red", pch=16)
#' curve(metadata(spk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # Not using spike-ins.
#' nspk <- modelGeneCV2(example.sce, use.spikes=FALSE)
#' nspk
#' 
#' plot(nspk$mean, nspk$cv2, pch=16, log="xy")
#' curve(metadata(nspk)$trend(x), add=TRUE, col="dodgerblue")
#'
#' # With blocking (and spike-ins).
#' block <- sample(LETTERS[1:2], ncol(example.sce), replace=TRUE)
#' blk <- modelGeneCV2(example.sce, block=block)
#' blk
#'
#' par(mfrow=c(1,2))
#' for (i in colnames(blk$per.block)) {
#'     current <- blk$per.block[[i]]
#'     plot(current$mean, current$cv2, pch=16, log="xy")
#'     points(metadata(current)$mean, metadata(current)$cv2, col="red", pch=16)
#'     curve(metadata(current)$trend(x), add=TRUE, col="dodgerblue")
#' }
#' 
#' @name modelGeneCV2
#' @aliases modelGeneCV2 modelGeneCV2,ANY-method modelGeneCV2,SingleCellExperiment-method
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm 
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom Matrix rowMeans
.model_gene_cv2_per_block <- function(x, ..., fit.size.factors=size.factors, subset.fit=NULL, 
    size.factors=NULL, subset.row=NULL) 
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x)
    }
    fit <- fitTrendCV2(x, size.factors=fit.size.factors, ..., subset.row=subset.fit)

    if (identical(subset.fit, subset.row)) {
        x.mean <- fit$mean
        x.cv2 <- fit$cv2
    } else {
        subset.row <- .subset_to_index(subset.row, x)
        x <- sweep(x[subset.row,,drop=FALSE], 2, size.factors, "/")
        x.mean <- rowMeans(x)
        x.cv2 <- rowVars(x)/x.mean^2
    }

    output <- DataFrame(mean=x.mean, cv2=x.cv2, trend=fit$trend(x.mean))
    output$ratio <- output$cv2/output$trend
    output$p.value <- pnorm(output$ratio, mean=1, sd=fit$std.dev, lower.tail=FALSE)

    rownames(output) <- rownames(x)[.subset_to_index(subset.row, x)]
    metadata(output) <- fit

    output
}

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pnorm p.adjust
.model_gene_cv2 <- function(x, ..., fit.size.factors=NULL, subset.fit=NULL, size.factors=NULL, subset.row=NULL, 
    block=NULL, equiweight=TRUE, method="fisher")
{
    args <- list(..., subset.fit=subset.fit, subset.row=subset.row)

    if (is.null(block)) {
        xargs <- list(x=x, fit.size.factors=fit.size.factors, size.factors=size.factors)
        output <- do.call(.model_gene_cv2_per_block, c(xargs, args))
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
        cur.args <- list(
            x=x[,current], 
            fit.size.factors=fit.size.factors[current],
            size.factors=size.factors[current]
        )
        collected[[i]] <- do.call(.model_gene_cv2_per_block, c(cur.args, args))
    }

    # Combining statistics with optional weighting.
    # We taking the geometric mean, which is more natural when the final statistic is a ratio.
    combined <- list()
    for (i in c("mean", "cv2", "trend", "ratio")) {
        extracted <- lapply(collected, "[[", i=i)
        extracted <- lapply(extracted, log)
        extracted <- mapply("*", extracted, weights, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        combined[[i]] <- exp(Reduce("+", extracted)/sum(weights))
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
#' @rdname modelGeneCV2
setGeneric("modelGeneCV2", function(x, ...) standardGeneric("modelGeneCV2"))

#' @export
#' @rdname modelGeneCV2
setMethod("modelGeneCV2", "ANY", .model_gene_cv2)

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment isSpike
#' @rdname modelGeneCV2
setMethod("modelGeneCV2", "SingleCellExperiment", function(x, ..., subset.fit=NULL, assay.type="counts", use.spikes=TRUE) {
    size.factors <- sizeFactors(x)
    if (use.spikes) {
        subset.fit <- isSpike(x)
        if (!any(subset.fit)) {
            stop("no spike-ins to use with 'use.spikes=TRUE'")
        }
        fit.size.factors <- sizeFactors(x, spikeNames(x)[1])
    } else {
        fit.size.factors <- size.factors
    }

    .check_centered_SF(x, assay.type=assay.type)
    .model_gene_cv2(assay(x, i=assay.type), ..., subset.fit=subset.fit, 
        fit.size.factors=fit.size.factors, size.factors=size.factors)
}) 
