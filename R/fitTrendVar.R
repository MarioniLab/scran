#' Fit a trend to the variances of log-counts
#' 
#' Fit a mean-dependent trend to the variances of the log-normalized expression values derived from count data.
#'
#' @param means A numeric vector containing the mean log-expression value for each gene.
#' @param vars A numeric vector containing the variance of log-expression values for each gene.
#' @param parametric A logical scalar indicating whether a parametric fit should be attempted.
#' @param lowess A logical scalar indicating whether a LOWESS fit should be attempted.
#' @param density.weights A logical scalar indicating whether to use inverse density weights.
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
#' The default fitting process follows a two-step procedure when \code{parametric=TRUE} and \code{lowess=TRUE}:
#' \enumerate{
#' \item A non-linear curve of the form
#' \deqn{y = \frac{ax}{x^n + b}}{y = ax/(x^n + b)}
#' is fitted to the variances against the means using \code{\link{nls}}.
#' Starting values and the number of iterations are automatically set if not explicitly specified in \code{nls.args}.
#' \item \code{\link{weightedLowess}} is applied to the log-ratios of the variance of each gene to the corresponding fitted value from the non-linear curve.
#' The final trend is defined as the product of the fitted values from the non-linear curve and the exponential of the LOWESS fitted value. 
#' If any tuning is necessary, the most important parameter here is \code{span}, which can be passed in the \code{...} arguments.
#' }
#' The aim is to use the parametric curve to reduce the sharpness of the expected mean-variance relationship for easier smoothing.
#' Conversely, the parametric form is not exact, so the smoothers will model any remaining trends in the residuals.
#'
#' If \code{parametric=FALSE}, LOWESS smoothing is performed directly on the log-variances using \code{\link{weightedLowess}}.
#' This may be helpful in situations where the data does not follow the expected parametric curve,
#' e.g., for transformed data that spans negative values where the expression is not defined.
#' (Indeed, the expression above is purely empirical, chosen simply as it matched the shape of the observed trend in many datasets.)
#'
#' If \code{lowess=FALSE}, the LOWESS smoothing step is omitted and the parametric fit is used directly.
#' This may be necessary in situations where the LOWESS overfits, e.g., because of very few points at high abundances.
#'
#' @section Filtering by mean:
#' Genes with mean log-expression below \code{min.mean} are not used in trend fitting.
#' This aims to remove the majority of low-abundance genes and preserve the sensitivity of span-based smoothers at moderate-to-high abundances.
#' It also protects against discreteness, which can interfere with estimation of the variability of the variance estimates and accurate scaling of the trend.
#'
#' Filtering is applied on the mean log-expression to avoid introducing spurious trends at the filter boundary.
#' The default threshold is chosen based on the point at which discreteness is observed in variance estimates from Poisson-distributed counts.
#' For heterogeneous droplet data, a lower threshold of 0.001-0.01 may be more appropriate, though this usually does not matter all too much.
#'
#' When extrapolating to values below the smallest observed mean (or \code{min.mean}), the output function will approach zero as the mean approaches zero.
#' This reflects the fact that the variance should be zero at a log-expression of zero (assuming a pseudo-count of 1 was used).
#' When extrapolating to values above the largest observed mean, the output function will be set to the fitted value of the trend at the largest mean.
#'
#' @section Weighting by density:
#' All fitting (with \code{\link{nls}} and \code{\link{weightedLowess}}) is performed by weighting each observation according to the inverse of the density of observations at the same mean.
#' This ensures that parts of the curve with few points are fitted as well as parts of the trend with many points.
#' Otherwise, differences in the distribution of means would favor good fits in highly dense intervals at the expense of sparser intervals.
#' (Note that these densities are computed after filtering on \code{min.mean}, so the high density of points at zero has no effect.)
#'
#' That being said, the density weighting can give inappropriate weight to very sparse intervals, especially those at high abundance.
#' This results in overfitting where the trend is compelled to pass through each point at these intervals.
#' For most part, this is harmless as most high-abundance genes are not highly variable so an overfitted trend is actually appropriate.
#' However, if high-abundance genes are variable, it may be better to set \code{density.weights=FALSE} to avoid this overfitting effect.
#' 
#' @author Aaron Lun
#' @seealso
#' \code{\link{modelGeneVar}} and \code{\link{modelGeneVarWithSpikes}}, where this function is used.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Fitting a trend:
#' library(DelayedMatrixStats)
#' means <- rowMeans(logcounts(sce))
#' vars <- rowVars(logcounts(sce))
#' fit <- fitTrendVar(means, vars)
#' 
#' # Comparing the two trend fits:
#' plot(means, vars, pch=16, cex=0.5, xlab="Mean", ylab="Variance")
#' curve(fit$trend(x), add=TRUE, col="dodgerblue", lwd=3)
#'
#' @export
#' @importFrom limma weightedLowess
#' @importFrom stats nls predict fitted approxfun
fitTrendVar <- function(means, vars, min.mean=0.1, parametric=TRUE, lowess=TRUE, density.weights=TRUE, nls.args=list(), ...) {
    # Filtering out zero-variance and low-abundance genes.
    is.okay <- !is.na(vars) & vars > 1e-8 & means >= min.mean 
    v <- vars[is.okay]
    m <- means[is.okay]

    if (density.weights) {
        w <- .inverse_density_weights(m, adjust=1)
    } else {
        w <- rep(1, length(m))
    }

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

        nls.args <- .setup_nls_args(nls.args, start.args=list(vars=v, means=m))
        nls.args$formula <- v ~ (exp(A)*m)/(m^(1+exp(N)) + exp(B))
        nls.args$weights <- w
        nls.args$control$warnOnly <- FALSE

        init.fit <- try(do.call(nls, nls.args), silent=TRUE) 
        if (is(init.fit, "try-error")) {
            Aest <- exp(nls.args$start$A)
            Best <- exp(nls.args$start$B)
            Nest <- exp(nls.args$start$N)+1
            PARAMFUN <- function(x) { Aest * x / (x^Nest + Best) }
            to.fit <- to.fit - log(PARAMFUN(m))
        } else {
            to.fit <- to.fit - log(fitted(init.fit))
            PARAMFUN <- function(x) { predict(init.fit, data.frame(m=x)) }
        }
    } else if (!lowess) {
        stop("at least one of 'lowess' or 'parametric' must be 'TRUE'")
    }

    if (lowess) {
        lfit <- weightedLowess(m, to.fit, weights=w, ...)
        LOESSFUN <- approxfun(m, lfit$fitted, rule=2)

        # Only trusting the parametric function for extrapolation; 
        # restricting non-parametric forms within the supported range.
        UNSCALEDFUN <- function(x) { 
            exp(LOESSFUN(x)) * PARAMFUN(x)
        }
    } else {
        UNSCALEDFUN <- PARAMFUN
    }

    # Adjusting for any scale shift due to fitting to the log-values.
    .correct_logged_expectation(m, v, w, UNSCALEDFUN)
}

#########################################################
# Computing NLS starting points for parametric fitting. #
#########################################################

#' @importFrom stats nls nls.control
.setup_nls_args <- function(nls.args, start.args) {
    if (is.null(nls.args)) {
        return(list())
    }

    nls.call <- do.call(call, c("nls", nls.args))
    nls.call <- match.call(nls, nls.call)
    nls.args <- as.list(nls.call)[-1]

    control <- nls.control(warnOnly=TRUE, maxiter=500)
    raw.start <- do.call(.get_nls_starts, start.args)
    start <- list(A=log(raw.start$a),
                  B=log(raw.start$b),
                  N=log(pmax(1e-8, raw.start$n-1))) # reflects enforced positivity in formula.

    altogether <- c(nls.args, list(control=control, start=start))
    keep <- !duplicated(names(altogether)) # nls.args are favoured.
    return(altogether[keep])
}

#' @importFrom stats coef lm fitted
#' @importFrom utils head tail
.get_nls_starts <- function(vars, means, left.n=100, left.prop=0.1,
    grid.length=10, b.grid.range=5, n.grid.max=10)
{
    o <- order(means)
    n <- length(vars)

    # Estimating the gradient from the left.
    left.n <- min(left.n, n*left.prop)
    keep <- head(o, max(1, left.n))
    y <- vars[keep]
    x <- means[keep]
    grad <- coef(lm(y~0+x))

    # Two-dimensional grid search is the most reliable way of estimating the remaining parameters.
    b.grid.pts <- 2^seq(from=-b.grid.range, to=b.grid.range, length.out=grid.length)
    n.grid.pts <- 2^seq(from=0, to=n.grid.max, length.out=grid.length)
    hits <- expand.grid(B=b.grid.pts, n=n.grid.pts)

    grid.ss <- mapply(B=hits$B, n=hits$n, FUN=function(B, n) {
        resid <- vars - (grad*B*means)/(means^n + B)
        sum(resid^2)
    })

    chosen <- which.min(grid.ss)
    N <- hits$n[chosen]
    B <- hits$B[chosen]
    A <- B * grad
    list(n=N, b=B, a=A)
}


