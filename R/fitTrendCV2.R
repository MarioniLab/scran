#' Fit a trend to the CV2
#' 
#' Fit a mean-dependent trend to the squared coefficient of variation, 
#' computed from count data after size factor normalization.
#'
#' @param means A numeric vector containing mean normalized expression values for all genes.
#' @param cv2 A numeric vector containing the squared coefficient of variation computed from normalized expression values for all genes.
#' @param ncells Integer scalar specifying the number of cells used to compute \code{cv2} and \code{means}.
#' @param min.mean Numeric scalar specifying the minimum mean to use for trend fitting.
#' @param top.prop A numeric scalar in (0, 1), specifying the top percentage of high-abundance genes to use for estimating the initial value of the plateau parameter of the fitted trend.
#' @param nls.args A list of parameters to pass to \code{\link{nls}}.
#' @param max.iter Integer scalar specifying the maximum number of robustness iterations to perform.
#' @param nmads Numeric scalar specifying the number of MADs to use to compute the tricube bandwidth during robustification.
#'
#' @return 
#' A named list is returned containing:
#' \describe{
#' \item{\code{trend}:}{A function that returns the fitted value of the trend at any value of the mean.}
#' \item{\code{std.dev}:}{A numeric scalar containing the robust standard deviation of the ratio of \code{var} to the fitted value of the trend across all features used for trend fitting.}
#' }
#'
#' @details
#' This function fits a mean-dependent trend to the CV2 of normalized expression values for the selected features.
#' Specifically, it fits a trend of the form
#' \deqn{y = A + \frac{B}{x}}{y = A + B/x}
#' using \code{\link{nls}}.
#' An initial value for \eqn{A} is estimated by obtaining the location of the asymptote at large means,
#' using \code{top.prop} to obtain a definition of \dQuote{large}.
#'
#' Trend fitting is performed after weighting each observation according to the inverse of the density of observations at the same mean.
#' This avoids problems with differences in the distribution of means that would otherwise favor good fits in highly dense intervals at the expense of sparser intervals.
#' Low-abundance genes with mean log-expression below \code{min.mean} are also removed prior to fitting, to avoid problems with discreteness and the upper bound on the CV2 at low counts.
#' 
#' Robustness iterations are also performed to protect against outliers. 
#' An initial fit is performed and each observation is weighted using tricube-transformed residuals (in addition to the existing inverse-density weights).
#' The bandwidth of the tricube scheme is defined as \code{nmads} multiplied by the median residual.
#' Iterations are performed until convergence or \code{max.iters} is reached.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(scater)
#' data(example.sce)
#' normcounts <- normalizeCounts(example.sce, log=FALSE)
#'
#' # Fitting a trend:
#' library(DelayedMatrixStats)
#' means <- rowMeans(normcounts)
#' cv2 <- rowVars(normcounts)/means^2
#' fit <- fitTrendCV2(means, cv2, ncol(example.sce))
#' 
#' # Examining the trend fit:
#' plot(means, cv2, pch=16, cex=0.5,
#'     xlab="Mean", ylab="CV2", log="xy")
#' curve(fit$trend(x), add=TRUE, col="dodgerblue", lwd=3)
#'
#' @seealso
#' \code{\link{modelGeneCV2}} and \code{\link{modelGeneCV2WithSpikes}}, where this function is used.
#' @export
#' @importFrom stats nls median quantile coef
fitTrendCV2 <- function(means, cv2, ncells, min.mean=0.1, top.prop=0.01, nls.args=list(), nmads=6, max.iter=50)
# Fits a spline to the log-CV2 values.
#
# written by Aaron Lun
# created 9 February 2017
{
    # Ignoring maxed CV2 values due to an outlier (caps at the number of cells).
    # Also ignoring means that are too low.
    to.use <- cv2 < ncells - 1e-8 & means > min.mean

    y <- cv2[to.use]
    x <- means[to.use]
    ly <- log(y)
    w <- .inverse_density_weights(log(x))

    # Rough estimation of initial parameters.
    # This is based on what happens as X->Inf (for A) or X->0 (for B).
    B <- median(ly + log(x), na.rm=TRUE)
    A <- median(ly[x >= quantile(x, 1-top.prop)])

    predFUN <- function(fit) {
        Aest <- exp(coef(fit)["A"])
        Best <- exp(coef(fit)["B"])
        function(x) { Aest  + Best/x }
    }

    # Robustness iterations until convergence. 
    # Weights are computed using residuals in the *log* space,
    # as we are currently fitting on the log-CV2 values.
    weights <- w
    tol <- 1e-8

    nls.args$formula <- ly ~ log(exp(A) + exp(B)/x)
    nls.args$start <- list(A=A, B=B)

    for (i in seq_len(max.iter)) {
        nls.args$weights <- weights
        fit <- do.call(nls, nls.args)

        r <- abs(ly - log(predFUN(fit)(x)))
        r <- r/(median(r, na.rm=TRUE) * nmads)
        r <- pmin(r, 1)

        new.weights <- (1 - r^3)^3 * w
        if (max(abs(new.weights - weights)) < tol) {
            break
        }
        weights <- new.weights
    }

    .correct_logged_expectation(x, y, w, predFUN(fit))
}
