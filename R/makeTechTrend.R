#' @importFrom BiocGenerics sizeFactors counts
#' @importFrom Matrix rowMeans
#' @importFrom SingleCellExperiment logcounts 
#' @importFrom BiocGenerics sizeFactors
#' @importFrom stats splinefun
#' @importFrom scater librarySizeFactors
#' @importFrom BiocParallel SerialParam bplapply bpmapply
#' @export
makeTechTrend <- function(means, size.factors=1, tol=1e-6, dispersion=0, pseudo.count=1, approx.npts=Inf, x=NULL, BPPARAM=SerialParam()) 
# This generates NB-distributed counts with the specified 
# dispersion in order to fit the mean-variance trend to the
# log-normalized counts. Designed for droplet data where
# spike-ins cannot be added to measure technical noise.
#
# written by Aaron Lun
# created 2 January 2018
{
    if (!is.null(x)) {
        .check_centered_SF(x, "logcounts")

        size.factors <- sizeFactors(x)
        if (is.null(size.factors)) { 
            size.factors <- librarySizeFactors(x, exprs_values="counts")
        }

        pseudo.count <- .get_log_offset(x)
        if (is.null(pseudo.count)) {
            stop("'log.exprs.offset' not specified in 'metadata(x)'")
        }

        all.ave <- rowMeans(logcounts(x))
        upper.value <- max(all.ave)
        means <- 2^seq(from=0, to=upper.value, length.out=100) - pseudo.count
    }

    to.core <- .worker_assign(length(means), BPPARAM)
    by.core <- .split_vector_by_workers(means, to.core)

    if (is.finite(approx.npts)) {
        if (approx.npts < 2) {
            stop("'approx.npts' should be at least 2")
        }

        # Approximating by fitting a spline with respect to the size factors.
        # This avoids having to compute these statistics for each size factor.
        lsf <- log(size.factors)
        if (length(lsf) <= approx.npts) {
            lpts <- lsf
        } else {
            lpts <- seq(min(lsf), max(lsf), length.out=approx.npts)
        }
        pts <- exp(lpts)

        out.expected <- bplapply(by.core, FUN=.tech_mean_computer, size.factors=pts, 
            tol=tol, dispersion=dispersion, pseudo.count=pseudo.count, BPPARAM=BPPARAM)
        collected.means <- .interpolate_and_average(lpts, unlist(out.expected, recursive=FALSE), lsf)

        by.core.constant <- .split_vector_by_workers(collected.means, to.core)
        out.sqdiff <- bpmapply(FUN=.tech_var_computer, means=by.core, logmeans=by.core.constant, 
            MoreArgs=list(size.factors=pts, tol=tol, dispersion=dispersion, pseudo.count=pseudo.count), 
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        collected.vars <- .interpolate_and_average(lpts, unlist(out.sqdiff, recursive=FALSE), lsf)

    } else {
        # Calling the C++ code to do the heavy lifting.
        raw.out <- bplapply(by.core, FUN=.tech_trend_computer, size.factors=size.factors, 
            tol=tol, dispersion=dispersion, pseudo.count=pseudo.count, BPPARAM=BPPARAM)

        collected.means <- unlist(lapply(raw.out, "[[", i=1))
        collected.vars <- unlist(lapply(raw.out, "[[", i=2))
    }

    # Creating a spline output function.
    splinefun(collected.means, collected.vars)
}

.tech_trend_computer <- function(means, size.factors, tol, dispersion, pseudo.count) 
# A helper function to ensure that the scran namespace is used in SnowParam().
{
    .Call(cxx_calc_log_count_stats, means, size.factors, tol, dispersion, pseudo.count)
}

.tech_mean_computer <- function(means, size.factors, tol, dispersion, pseudo.count)
# A helper function to ensure that the scran namespace is used in SnowParam().
{
    .Call(cxx_calc_log_expected, means, size.factors, tol, dispersion, pseudo.count)
}

.tech_var_computer <- function(means, size.factors, tol, dispersion, pseudo.count, logmeans)
# A helper function to ensure that the scran namespace is used in SnowParam().
{
    .Call(cxx_calc_log_sqdiff, means, size.factors, tol, dispersion, pseudo.count, logmeans)
}

#' @importFrom stats spline
.interpolate_and_average <- function(x, y.list, xout) 
# Fits a spline to x against each element of y,
# interpolates at xout and then computes the average.
{
    vapply(y.list, FUN=function(y) {
        mean(spline(x, y, xout=xout)$y)
    }, FUN.VALUE=0)
}
