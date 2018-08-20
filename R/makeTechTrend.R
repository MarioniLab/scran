#' @importFrom BiocGenerics sizeFactors counts
#' @importFrom Matrix colSums rowMeans
#' @importFrom SingleCellExperiment logcounts
#' @importFrom stats splinefun
#' @importFrom BiocParallel SerialParam bplapply
#' @export
makeTechTrend <- function(means, size.factors=1, tol=1e-6, dispersion=0, pseudo.count=1, x=NULL, BPPARAM=SerialParam()) 
# This generates NB-distributed counts with the specified 
# dispersion in order to fit the mean-variance trend to the
# log-normalized counts. Designed for droplet data where
# spike-ins cannot be added to measure technical noise.
#
# written by Aaron Lun
# created 2 January 2018
{
    if (!is.null(x)) {
        size.factors <- sizeFactors(x)
        if (is.null(size.factors)) { 
            size.factors <- colSums(counts(x))
            size.factors <- size.factors/mean(size.factors)
        }

        pseudo.count <- .get_log_offset(x)
        if (is.null(pseudo.count)) {
            stop("'log.exprs.offset' not specified in 'metadata(x)'")
        }

        all.ave <- rowMeans(logcounts(x))
        upper.value <- max(all.ave)
        means <- 2^seq(from=0, to=upper.value, length.out=100) - pseudo.count
    }

    if (abs(mean(size.factors) - 1) > 1e-6) {
        stop("size factors should be centred at unity") 
    }

    # Calling the C++ code to do the heavy lifting.
    to.core <- .worker_assign(length(means), BPPARAM)
    by.core <- .split_vector_by_workers(means, to.core)
    raw.out <- bplapply(by.core, FUN=.tech_trend_computer, size.factors=size.factors, 
        tol=tol, dispersion=dispersion, pseudo.count=pseudo.count, BPPARAM=BPPARAM)

    collected.means <- unlist(lapply(raw.out, "[[", i=1))
    collected.vars<- unlist(lapply(raw.out, "[[", i=2))

    # Creating a spline output function.
    splinefun(collected.means, collected.vars)
}

.tech_trend_computer <- function(means, size.factors, tol, dispersion, pseudo.count) 
# A helper function to ensure that the scran namespace is used in SnowParam().
{
    .Call(cxx_calc_log_count_stats, means, size.factors, tol, dispersion, pseudo.count)
}
