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
    args <- list(tol=tol, disp=disp, pseudo=pseudo)

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
        args$Sizes <- pts

        out.expected <- bpmapply(Means=by.core, FUN=calc_log_expected, MoreArgs=args,
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        collected.means <- .interpolate_and_average(lpts, unlist(out.expected, recursive=FALSE), lsf)

        by.core.constant <- .split_vector_by_workers(collected.means, to.core)
        out.sqdiff <- bpmapply(FUN=calc_log_sqdiff, Means=by.core, Constants=by.core.constant, MoreArgs=args,
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        collected.vars <- .interpolate_and_average(lpts, unlist(out.sqdiff, recursive=FALSE), lsf)

    } else {
        args$Sizes <- size.factors
        raw.out <- bpmapply(Means=by.core, FUN=calc_log_count_stats, MoreArgs=args,
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

        collected.means <- unlist(lapply(raw.out, "[[", i=1))
        collected.vars <- unlist(lapply(raw.out, "[[", i=2))
    }

    # Creating a spline output function.
    splinefun(collected.means, collected.vars)
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
