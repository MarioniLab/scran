.ranksafe_qr <- function(design, tol=1e-7) 
# Rank-checking QR decomposition of a design matrix. Throws an
# error if the design matrix is not of full rank, which simplifies
# downstream processes as full rank can always be assumed.
{
    out <- qr(design, LAPACK=TRUE)
    d <- diag(out$qr)
    if (!all(abs(d) > tol)) { 
        stop("design matrix is not of full rank")
    }
    return(out)
}

.calc_residuals_wt_zeroes <- function(x, design, QR, subset.row, lower.bound) 
# Computes residuals, but ensures that residuals for observations that 
# were below the lower bound are set to a constant value that is smaller 
# than all other residuals. Avoids spurious patterns when
# 
{
    if (!missing(design)) {
        QR <- .ranksafe_qr(design)
    }
    if (is.null(lower.bound)) { 
        stop("lower bound must be supplied or NA when computing residuals")
    }
    get_residuals(x, QR$qr, QR$qraux, subset.row - 1L, as.double(lower.bound))
}

.guess_lower_bound <- function(x, assay.type, lower.bound) 
# Getting the lower bound on the expression values for a given assay, if not supplied.
# We bump it up a little to make sure that expression values at the lower bound will 
# actually be detected as being "<= lower.bound".
{ 
    if (is.null(lower.bound)) { 
        if (assay.type=="logcounts") {
            lower.bound <- log2(.get_log_offset(x)) + 1e-8
        } else if (assay.type=="counts") {
            lower.bound <- 1e-8
        }
    }
    return(lower.bound)
}

#' @importFrom S4Vectors metadata
.get_log_offset <- function(x) 
# Helper function to get the log-offset value.
{
    metadata(x)$log.exprs.offset
}
