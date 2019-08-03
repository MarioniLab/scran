#' Fit a trend to the variances of log-counts
#' 
#' Fit a mean-dependent trend to the variances of the log-normalized expression values derived from count data.
#'
#' @param x A numeric matrix of log-counts, or a \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param parametric A logical scalar indicating whether a parametric fit should be attempted.
#' @param nls.args A list of parameters to pass to \code{\link{nls}} if \code{parametric=TRUE}.
#' @param design A numeric matrix containing blocking terms for uninteresting factors of variation.
#' @param min.mean A numeric scalar specifying the minimum mean to use for trend fitting.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param ... For the \code{ANY} method, further arguments to pass to \code{\link{weightedLowess}} for LOWESS fitting.
#'
#' For the generic and SingleCellExperiment methods, further arguments to pass to the \code{ANY} method.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether parallelization should be performed across genes.
#'
#' @return 
#' A named list is returned, containing:
#' \describe{
#' \item{\code{mean}:}{A numeric vector of mean log-expression values for all features used for trend fitting.}
#' \item{\code{var}:}{A numeric vector of the variances of log-expression values for all features used for trend fitting.}
#' \item{\code{trend}:}{A function that returns the fitted value of the trend at any value of the mean.}
#' \item{\code{std.dev}:}{A numeric scalar containing the robust standard deviation of the ratio of \code{var} to the fitted value of the trend across all features used for trend fitting.}
#' }
#'
#' @details
#' This function fits a mean-dependent trend to the variance of the log-normalized expression for the selected features.
#' The fitted trend can then be used to decompose the variance of each gene into biological and technical components.
#' Log-transformed values are used as these are more robust to genes/transcripts with strong expression in a few outlier cells.
#'
#' @author Aaron Lun
#'
#' @examples
#' data(example.sce)
#'
#' # Fitting a trend:
#' fit.spike <- fitTrendVar(example.sce) # defaults to using spike-ins
#' fit.endog <- fitTrendVar(example.sce, use.spikes=FALSE)
#' 
#' # Comparing the two trend fits:
#' plot(fit.endog$mean, fit.endog$var, pch=16, cex=0.5,
#'     xlab="Mean", ylab="Variance")
#' points(fit.spike$mean, fit.spike$var, col="red", pch=16)
#' 
#' curve(fit.endog$trend(x), add=TRUE, col="dodgerblue", lwd=3)
#' curve(fit.spike$trend(x), add=TRUE, col="salmon", lwd=3)
#'
#' @name fitTrendVar
#' @aliases fitTrendVar fitTrendVar,ANY-method fitTrendVar,SingleCellExperiment-method
NULL

#############################
# Defining the basic method #
#############################

#' @importFrom BiocParallel SerialParam
.fit_trend_var <- function(x, parametric=TRUE, nls.args=list(), ..., design=NULL, 
    min.mean=0.1, subset.row=NULL, BPPARAM=SerialParam())
{
    stats.out <- .lognormvar(x, block=NULL, design=design, subset.row=subset.row, BPPARAM=BPPARAM)
    .fit_trend_var0(stats.out$means, stats.out$vars, min.mean=min.mean, parametric=parametric, nls.args=nls.args, ...)
}

#' @importFrom stats nls fitted median predict
#' @importFrom limma weightedLowess weighted.median
.fit_trend_var0 <- function(means, vars, min.mean=0.1, parametric=TRUE, nls.args=list(), ...) {
    # Filtering out zero-variance and low-abundance genes.
    is.okay <- !is.na(vars) & vars > 1e-8 & means >= min.mean 
    v <- vars[is.okay]
    m <- means[is.okay]
    w <- .inverse_density_weights(m, adjust=1)

    # Default parametric trend is a straight line from 0 to 1 below the supported range.
    if (length(v) < 2L) {
        stop("need at least 2 points for non-parametric curve fitting")
    } 
    to.fit <- log(v)
    left.edge <- min(m)
    PARAMFUN <- function(x) { pmin(1, x/left.edge) } 

    # Fitting a parametric curve to try to flatten the shape.
    # This is of the form y = a*x/(x^n + b), but each coefficent is actually set
    # to exp(*) to avoid needing to set lower bounds.
    if (parametric) { 
        attempt <- try({
            if (length(v) <= 3L) {
                stop("need at least 4 points for non-linear curve fitting")
            } 

            nls.args <- .setup_nls_args(nls.args, start.args=list(vars=v, means=m))
            nls.args$formula <- v ~ (exp(A)*m)/(m^(1+exp(N)) + exp(B))
            nls.args$weights <- w

            init.fit <- do.call(nls, nls.args)
            to.fit <- to.fit - log(fitted(init.fit))
            PARAMFUN <- function(x) { predict(init.fit, data.frame(m=x)) }
        }, silent=TRUE)

        if (is(attempt, "try-error")) {
            warning("parametric curve fitting failed, defaulting to loess-only")
        }
    } 

    lfit <- weightedLowess(m, to.fit, weights=w, ...)
    LOESSFUN <- approxfun(m, lfit$fitted, rule=2)

    # Only trusting the parametric function for extrapolation; 
    # restricting non-parametric forms within the supported range.
    UNSCALEDFUN <- function(x) { 
        exp(LOESSFUN(x)) * PARAMFUN(x)
    }

    # Adjusting for any scale shift due to fitting to the log-values.
    corrected <- .correct_logged_expectation(m, v, w, UNSCALEDFUN)

    c(list(mean=means, var=vars), corrected)
}

#########################
# Setting up S4 methods #
#########################

#' @export
#' @rdname fitTrendVar
setGeneric("fitTrendVar", function(x, ...) standardGeneric("fitTrendVar"))

#' @export
#' @rdname fitTrendVar
setMethod("fitTrendVar", "ANY", .fit_trend_var)

#' @export
#' @rdname fitTrendVar
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment isSpike
setMethod("fitTrendVar", "SingleCellExperiment", function(x, subset.row=NULL, ..., assay.type="logcounts", use.spikes=TRUE) {
    .check_centered_SF(x, assay.type=assay.type)
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Can use only spikes, everything but spikes, or everything. 
    if (!is.na(use.spikes)) {
        is.spike <- isSpike(x)
        if (is.null(is.spike)) {
            is.spike <- logical(nrow(x))
        }
        is.spike <- which(is.spike)
        if (use.spikes) {
            if (!length(is.spike)) {
                warning("no spike-in transcripts present for 'use.spikes=TRUE'")
            } else {
                subset.row <- intersect(subset.row, is.spike)
            }
        } else {
            subset.row <- setdiff(subset.row, is.spike)
        }
    }

    .fit_trend_var(assay(x, i=assay.type), ..., subset.row=subset.row)
})
