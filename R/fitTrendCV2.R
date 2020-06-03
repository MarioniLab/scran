#' Fit a trend to the CV2
#' 
#' Fit a mean-dependent trend to the squared coefficient of variation, 
#' computed from count data after size factor normalization.
#'
#' @param means A numeric vector containing mean normalized expression values for all genes.
#' @param cv2 A numeric vector containing the squared coefficient of variation computed from normalized expression values for all genes.
#' @param ncells Integer scalar specifying the number of cells used to compute \code{cv2} and \code{means}.
#' @param min.mean Numeric scalar specifying the minimum mean to use for trend fitting.
#' @param nls.args A list of parameters to pass to \code{\link{nls}}.
#' @param max.iter Integer scalar specifying the maximum number of robustness iterations to perform.
#' @param nmads Numeric scalar specifying the number of MADs to use to compute the tricube bandwidth during robustification.
#' @param simplified Logical scalar indicating whether the function can automatically use a simpler trend if errors are encountered for the usual paramterization.
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
#' using an iteratively reweighted least-squares approach implemented via \code{\link{nls}}.
#' This trend is based on a similar formulation from \pkg{DESeq2} and generally captures the mean-CV2 trend well.
#'
#' Trend fitting is performed after weighting each observation according to the inverse of the density of observations at the same mean.
#' This avoids problems with differences in the distribution of means that would otherwise favor good fits in highly dense intervals at the expense of sparser intervals.
#' Low-abundance genes with means below \code{min.mean} are also removed prior to fitting, to avoid problems with discreteness and the upper bound on the CV2 at low counts.
#' 
#' Robustness iterations are also performed to protect against outliers. 
#' An initial fit is performed and each observation is weighted using tricube-transformed standardized residuals (in addition to the existing inverse-density weights).
#' The bandwidth of the tricube scheme is defined as \code{nmads} multiplied by the median standardized residual.
#' Iterations are performed until convergence or \code{max.iters} is reached.
#'
#' Occasionally, there are not enough high-abundance points to uniquely determine the \eqn{A} parameter.
#' In such cases, the function collapses back to fitting a simpler trend
#' \deqn{y = \frac{B}{x}}{y = B/x}
#' to avoid errors about singular gradients in \code{\link{nls}}.
#' If \code{simplified=FALSE}, this simplification is not allowed and the error is directly reported.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' normcounts <- normalizeCounts(sce, log=FALSE)
#'
#' # Fitting a trend:
#' library(DelayedMatrixStats)
#' means <- rowMeans(normcounts)
#' cv2 <- rowVars(normcounts)/means^2
#' fit <- fitTrendCV2(means, cv2, ncol(sce))
#' 
#' # Examining the trend fit:
#' plot(means, cv2, pch=16, cex=0.5,
#'     xlab="Mean", ylab="CV2", log="xy")
#' curve(fit$trend(x), add=TRUE, col="dodgerblue", lwd=3)
#'
#' @references
#' Brennecke P, Anders S, Kim JK et al. (2013).
#' Accounting for technical noise in single-cell RNA-seq experiments.
#' \emph{Nat. Methods} 10:1093-95
#' 
#' @seealso
#' \code{\link{modelGeneCV2}} and \code{\link{modelGeneCV2WithSpikes}}, where this function is used.
#' @export
#' @importFrom stats nls median coef
#' @importFrom statmod glmgam.fit
fitTrendCV2 <- function(means, cv2, ncells, min.mean=0.1, nls.args=list(), 
    simplified=TRUE, nmads=6, max.iter=50)
{
    # Ignoring maxed CV2 values due to an outlier (caps at the number of cells).
    # Also ignoring means that are too low.
    to.use <- cv2 < ncells - 1e-8 & means > min.mean

    y <- cv2[to.use]
    x <- means[to.use]
    w <- .inverse_density_weights(log(x))

    # Rough estimation of initial parameters, assuming that var \propto mean^2
    # in terms of the distribution around the trend.
    out <- glmgam.fit(y=y, X=cbind(1, 1/x))
    coefs <- log(pmax(1e-8, out$coefficients))
    names(coefs) <- c("logA", "logB")

    predFUN <- function(coefs) {
        Aest <- unname(exp(coefs["logA"]))
        Best <- unname(exp(coefs["logB"]))
        function(x) Aest  + Best/x
    }

    # Robustness iterations until convergence. 
    # Every iteration should update 'coefs', 'fitted' and 'weights'.
    tol <- 1e-8
    weights <- w 
    fitted <- predFUN(coefs)(x)
    nls.args$formula <- y ~ exp(logA) + exp(logB)/x

    for (i in seq_len(max.iter)) {
        nls.args$start <- as.list(coefs)

        # Weights are always scaled by fitted values to adjust for mean-variance
        # relationship _of the estimates around the trend_. This means that we 
        # effectively perform IRLS via nls().
        nls.args$weights <- weights / fitted^2

        # nls() can fail to converge if the trend is just a straight line,
        # such that 'A' is any arbitrarily small value.
        fit <- try(do.call(nls, nls.args), silent=TRUE)
        if (is(fit, "try-error")) {
            msg <- attr(fit, "condition")$message
            if (simplified && grepl("singular gradient", msg)) {
                simplified <- FALSE # can't get here again.

                nls.args$formula <- y ~ exp(logB)/x
                nls.args$start$logA <- NULL
                fit <- do.call(nls, nls.args)
                predFUN <- function(coefs) {
                    Best <- unname(exp(coefs["logB"]))
                    function(x) Best/x
                }
            } else {
                stop(msg)
            }
        }

        coefs <- coef(fit)
        fitted <- predFUN(coefs)(x)
        r <- abs(y/fitted - 1)
        r <- r/(median(r, na.rm=TRUE) * nmads)
        r <- pmin(r, 1)

        new.weights <- (1 - r^3)^3 * w 
        if (max(abs(new.weights - weights)) < tol) {
            break
        }
        weights <- new.weights
    }

    FUN <- predFUN(coefs)
    leftovers <- y/FUN(x)
    std.dev <- unname(weighted.median(abs(leftovers - 1), w, na.rm=TRUE)) * 1.4826
    list(trend=FUN, std.dev=std.dev)
}
