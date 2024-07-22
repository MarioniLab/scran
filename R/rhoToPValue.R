#' Spearman's rho to a p-value
#'
#' Compute an approximate p-value against the null hypothesis that Spearman's rho is zero.
#' This vectorizes the approximate p-value calculation in \code{\link{cor.test}} with \code{method="spearman"}.
#'
#' @param rho Numeric vector of rho values.
#' @param n Integer scalar specifying the number of observations used to compute \code{rho}.
#' @param positive Logical scalar indicating whether to perform a one-sided test for the alternative of a positive (\code{TRUE}) or negative rho (\code{FALSE}).
#' Default is to return statistics for both directions.
#'
#' @return 
#' If \code{positive=NULL}, a list of two numeric vectors is returned,
#' containing p-values for the test against the alternative hypothesis in each direction.
#'
#' Otherwise, a numeric vector is returned containing the p-values for the test in the specified direction.
#'
#' @author Aaron Lun
#'
#' @examples
#' rhoToPValue(seq(-1, 1, 21), 50)
#'
#' @export
#' @importFrom stats pt
rhoToPValue <- function(rho, n, positive=NULL) {
    # Mildly adapted from cor.test.
    q <- (n^3 - n) * (1 - rho)/6
    den <- (n * (n^2 - 1)/6)
    r <- 1 - q/den
    tstat <- r/sqrt((1 - r^2)/(n - 2))

    FUN <- function(p) pt(tstat, df = n - 2, lower.tail = !p)

    if (!is.null(positive)) {
        FUN(positive)
    } else {
        list(positive=FUN(TRUE), negative=FUN(FALSE))
    }
}
