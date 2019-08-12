#' Generate a trend for Poisson noise
#' 
#' Create a mean-variance trend for log-normalized expression values derived from Poisson-distributed counts.
#'
#' @param means A numeric vector of length 2 or more, containing the range of mean counts observed in the dataset.
#' @param size.factors A numeric vector of size factors for all cells in the dataset.
#' @param dispersion A numeric scalar specifying the dispersion for the NB distribution.
#' If zero, a Poisson distribution is used.
#' @param pseudo.count A numeric scalar specifying the pseudo-count to be added to the scaled counts before log-transformation.
#' @param npts An integer scalar specifying the number of interpolation points to use.
#' @param ... Further arguments to pass to \code{\link{fitTrendVar}} for trend fitting.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how parallelization should be performed across interpolation points.
#'
#' @return A named list is returned containing:
#' \describe{
#' \item{\code{trend}:}{A function that returns the fitted value of the trend at any value of the mean.}
#' \item{\code{std.dev}:}{A numeric scalar containing the robust standard deviation of the ratio of \code{var} to the fitted value of the trend across all features used for trend fitting.}
#' }
#'
#' @details
#' This function is useful for modelling technical noise in highly diverse datasets without spike-ins,
#' where fitting a trend to the endogenous genes would not be appropriate given the strong biological heterogeneity.
#' It is mostly intended for UMI datasets where the technical noise is close to Poisson-distributed.
#' 
#' This function operates by simulating Poisson or negative binomial-distributed counts,
#' computing log-transformed normalized expression values from those counts,
#' calculating the mean and variance and then passing those metrics to \code{\link{fitTrendVar}}.
#' The log-transformation ensures that variance is modelled in the same space that is used for downstream analyses like PCA.
#'
#' Simulations are performed across a range of values in \code{means} to achieve reliable interpolation,
#' with the stability of the trend determined by the number of simulation points in \code{npts}.
#' The number of cells is determined from the length of \code{size.factors},
#' which are used to scale the distribution means prior to sampling counts.
#' 
#' @seealso
#' \code{\link{fitTrendVar}}, which is used to fit the trend.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Mocking up means and size factors:
#' sf <- 2^rnorm(1000, sd=0.1)
#' sf <- sf/mean(sf)
#' means <- rexp(100, 0.1)
#'
#' # Using these to construct a Poisson trend:
#' out <- fitTrendPoisson(means, sf)
#' curve(out$trend(x), xlim=c(0, 10))
#'
#' @export
#' @importFrom BiocParallel SerialParam
fitTrendPoisson <- function(means, size.factors, npts=1000, dispersion=0, pseudo.count=1, BPPARAM=SerialParam(), ...) {
    out <- .generate_poisson_values(means, size.factors, npts=npts, 
        dispersion=dispersion, pseudo.count=pseudo.count, BPPARAM=BPPARAM)
    fitTrendVar(out$means, out$vars, ...)
}

#' @importFrom stats rpois rnbinom
#' @importFrom BiocParallel SerialParam
.generate_poisson_values <- function(means, size.factors, npts=1000, dispersion=0, pseudo.count=1, 
    block=NULL, design=NULL, BPPARAM=SerialParam()) 
{
    if (dispersion==0) {
        FUN <- function(m) rpois(length(m), lambda=m)
    } else {
        FUN <- function(m) rnbinom(length(m), mu=m, size=1/dispersion)
    }

    means <- means[means>0]
    pts <- exp(seq(from=log(min(means)), to=log(max(means)), length=npts))

    Y <- matrix(0, npts, length(size.factors))
    for (i in seq_along(pts)) {
        Y[i,] <- FUN(pts[i] * size.factors)
    }

   .compute_mean_var(Y, block=block, design=design, subset.row=NULL,
        block.FUN=compute_blocked_stats_lognorm, 
        residual.FUN=compute_residual_stats_lognorm, 
        BPPARAM=BPPARAM, sf=size.factors, pseudo=pseudo.count)
}
