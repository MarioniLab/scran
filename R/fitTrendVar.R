#' Fit a trend to the variances of log-counts
#' 
#' Fit a mean-dependent trend to the variances of the log-normalized expression values derived from count data.
#'
#' @param means A numeric vector containing the mean log-expression value for each gene.
#' @param vars A numeric vector containing the variance of log-expression values for each gene.
#' @param parametric A logical scalar indicating whether a parametric fit should be attempted.
#' @param nls.args A list of parameters to pass to \code{\link{nls}} if \code{parametric=TRUE}.
#' @param min.mean A numeric scalar specifying the minimum mean to use for trend fitting.
#' @param ... Further arguments to pass to \code{\link{weightedLowess}} for LOWESS fitting.
#'
#' @return 
#' A named list is returned containing:
#' \describe{
#' \item{\code{trend}:}{A function that returns the fitted value of the trend at any value of the mean.}
#' \item{\code{std.dev}:}{A numeric scalar containing the robust standard deviation of the ratio of \code{var} to the fitted value of the trend across all features used for trend fitting.}
#' }
#'
#' @details
#' This function fits a mean-dependent trend to the variance of the log-normalized expression for the selected features.
#' The fitted trend can then be used to decompose the variance of each gene into biological and technical components,
#' as done in \code{\link{modelGeneVar}} and \code{\link{modelGeneVarWithSpikes}}.
#'
#' If \code{parametric=TRUE}, a non-linear curve of the form
#' \deqn{y = \frac{ax}{x^n + b}}{y = ax/(x^n + b)}
#' is fitted to the variances against the means using \code{\link{nls}}.
#' Starting values and the number of iterations are automatically set if not explicitly specified in \code{nls.args}.
#' \code{\link{weightedLowess}} is then applied to the log-ratios of the variance to the fitted value for each gene.
#' The aim is to use the parametric curve to reduce the sharpness of the expected mean-variance relationship[for easier smoothing.
#' Conversely, the parametric form is not exact, so the smoothers will model any remaining trends in the residuals.
#'
#' If \code{parametric=FALSE}, smoothing is performed directly on the log-variances using \code{\link{weightedLowess}}.
#'
#' Low-abundance genes with mean log-expression below \code{min.mean} are not used in trend fitting, to preserve the sensitivity of span-based smoothers at moderate-to-high abundances.
#' It also protects against discreteness, which can interfere with estimation of the variability of the variance estimates and accurate scaling of the trend.
#' The default threshold is chosen based on the point at which discreteness is observed in variance estimates from Poisson-distributed counts.
#' For heterogeneous droplet data, a lower threshold of 0.001-0.01 may be more appropriate.
#'
#' When extrapolating to values below the smallest observed mean (or \code{min.mean}), the output function will approach zero as the mean approaches zero.
#' This reflects the fact that the variance should be zero at a log-expression of zero (assuming a pseudo-count of 1 was used).
#' When extrapolating to values above the largest observed mean, the output function will be set to the fitted value of the trend at the largest mean.
#'
#' All fitting (with \code{\link{nls}} and \code{\link{weightedLowess}}) is performed by weighting each observation according to the inverse of the density of observations at the same mean.
#' This avoids problems with differences in the distribution of means that would otherwise favor good fits in highly dense intervals at the expense of sparser intervals.
#' Note that these densities are computed after filtering on \code{min.mean}.
#' 
#' @author Aaron Lun
#' @seealso
#' \code{\link{modelGeneVar}} and \code{\link{modelGeneVarWithSpikes}}, where this function is used.
#'
#' @examples
#' data(example.sce)
#'
#' # Fitting a trend:
#' library(DelayedMatrixStats)
#' means <- rowMeans(logcounts(example.sce))
#' vars <- rowVars(logcounts(example.sce))
#' fit <- fitTrendVar(means, vars)
#' 
#' # Comparing the two trend fits:
#' plot(means, vars, pch=16, cex=0.5, xlab="Mean", ylab="Variance")
#' curve(fit$trend(x), add=TRUE, col="dodgerblue", lwd=3)
#'
#' @export
#' @importFrom limma weightedLowess
#' @importFrom stats nls predict fitted approxfun
fitTrendVar <- function(means, vars, min.mean=0.1, parametric=TRUE, nls.args=list(), ...) {
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
        if (length(v) <= 3L) {
            stop("need at least 4 points for non-linear curve fitting")
        } 

        attempt <- try({
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
    .correct_logged_expectation(m, v, w, UNSCALEDFUN)
}
