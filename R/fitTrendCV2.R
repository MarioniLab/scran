#' Fit a trend to the CV2
#' 
#' Fit a mean-dependent trend to the squared coefficient of variation, 
#' computed from count data after size factor normalization.
#'
#' @param x A numeric matrix of counts, or a \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param size.factors A numeric vector of size factors for each cell/column in \code{x}.
#' @param top.prop A numeric scalar in (0, 1), specifying the top percentage of high-abundance genes to use for estimating the initial value of the plateau parameter of the fitted trend.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param ... For the \code{ANY} method, further arguments to pass to \code{\link{weightedLowess}} for LOWESS fitting.
#'
#' For the generic and SingleCellExperiment methods, further arguments to pass to the \code{ANY} method.
#'
#' @return 
#' A named list is returned, containing:
#' \describe{
#' \item{\code{mean}:}{A numeric vector of mean normalized expression values for all features used for trend fitting.}
#' \item{\code{cv2}:}{A numeric vector of the of CV2 values for all features used for trend fitting.}
#' \item{\code{trend}:}{A function that returns the fitted value of the trend at any value of the mean.}
#' \item{\code{std.dev}:}{A numeric scalar containing the robust standard deviation of the ratio of \code{var} to the fitted value of the trend across all features used for trend fitting.}
#' }
#'
#' @details
#' This function fits a mean-dependent trend to the CV2 of normalized expression values for the selected features.
#' The fitted trend can then be used to identify highly variable genes that have CV2 values above the trend.
#'
#' @author Aaron Lun
#'
#' @examples
#' data(example.sce)
#'
#' # Fitting a trend:
#' fit <- fitTrendCV2(example.sce) 
#' 
#' # Examining the trend fit:
#' plot(fit$mean, fit$cv2, pch=16, cex=0.5,
#'     xlab="Mean", ylab="CV2", log="xy")
#' curve(fit$trend(x), add=TRUE, col="dodgerblue", lwd=3)
#'
#' @export
#' @importFrom stats nls median quantile coef
fitTrendCV2 <- function(means, cv2, ncells, min.mean=0.1, subset.row=NULL, top.prop=0.01)
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
    B <- median(ly + log(x), na.rm=TRUE)
    A <- median(ly[x >= quantile(x, 1-top.prop)])

    fitFUN <- function(X, Y, W) {
        nls(Y ~ log(exp(A) + exp(B)/X), start=list(A=A, B=B), weights=W)
    }

    predFUN <- function(fit) {
        Aest <- exp(coef(fit)["A"])
        Best <- exp(coef(fit)["B"])
        function(x) { Aest  + Best/x }
    }

    # Robustness iterations until convergence. 
    # Weights are computed using residuals in the *log* space,
    # as we are currently fitting on the log-CV2 values.
    weights <- w
    nmads <- 6
    tol <- 1e-8
    max.iter <- 50

    for (i in seq_len(max.iter)) {
        fit <- fitFUN(x, ly, weights)

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
