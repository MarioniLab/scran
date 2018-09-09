#' @importFrom BiocGenerics sizeFactors counts
#' @importFrom Matrix colSums rowMeans
#' @importFrom SingleCellExperiment logcounts
#' @importFrom stats splinefun
#' @importFrom BiocParallel SerialParam bplapply
#' @export
makeTechTrend <- function(means, size.factors=1, tol=1e-6, dispersion=0, pseudo.count=1, approximate=FALSE, x=NULL, BPPARAM=SerialParam()) 
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

    if (approximate) {
        approx.out <- .approximate_log_stats(means, size.factors, dispersion, pseudo.count)
        collected.means <- approx.out$mean
        collected.vars <- approx.out$var

    } else {
        # Calling the C++ code to do the heavy lifting.
        to.core <- .worker_assign(length(means), BPPARAM)
        by.core <- .split_vector_by_workers(means, to.core)
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

.approximate_log_stats <- function(means, size.factors, dispersion, pseudo.count) 
# Applies Taylor series approximations to compute the moments.
{
    output.mean <- output.var <- numeric(length(means))
    for (i in seq_along(means)) {
        cur.mean <- means[i]
        cell.mean <- cur.mean * size.factors
        cell.var <- cell.mean + dispersion * cell.mean^2
        inside <- cur.mean + pseudo.count # equal to [E(X)/s + p]

        # E[log(X/s + p)] ~= log[E(X)/s + p] - 1/2 * var(X) / [E(X)/s + p]^2 / s^2 
        # ... where X ~ NB(mu*s, d)
        mus <- log(inside) - 0.5 * cell.var / size.factors^2 / inside^2
        grand.mean <- mean(mus)
        output.mean[i] <- grand.mean

        # E[log(X/s + p)^2] ~= log[E(X)/s + p]^2 + { 1/[E(X)/s + p]^2 - log[E(X)/s + p] / [E(X)/s + p]^2 } * var(X) / s^2
        exp.sq <- log(inside)^2 + ( 1/inside^2 - log(inside) / inside^2 ) * cell.var / size.factors^2

        # Remember that we want to compute the variance of the gene relative to the grand mean!
        output.var[i] <- mean(exp.sq - 2 * mus * grand.mean + grand.mean^2)
    }

    # Need to divide by log(2) or log(2)^2.
    output.mean <- output.mean/log(2)
    output.var <- output.var/log(2)^2
    list(mean=output.mean, var=output.var)
}


