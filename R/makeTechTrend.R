makeTechTrend <- function(means, size.factors=1, tol=1e-6, dispersion=0, pseudo.count=1, sce=NULL) 
# This generates NB-distributed counts with the specified 
# dispersion in order to fit the mean-variance trend to the
# log-normalized counts. Designed for droplet data where
# spike-ins cannot be added to measure technical noise.
#
# written by Aaron Lun
# created 2 January 2018
{
    if (!is.null(sce)) {
        size.factors <- sizeFactors(sce)
        if (is.null(size.factors)) { 
            size.factors <- colSums(counts(sce))
            size.factors <- size.factors/mean(size.factors)
        }

        pseudo.count <- .get_log_offset(sce)
        if (is.null(pseudo.count)) {
            stop("'log.exprs.offset' not specified in 'metadata(sce)'")
        }

        all.ave <- rowMeans(logcounts(sce))
        upper.value <- max(all.ave)
        means <- 2^seq(from=0, to=upper.value, length.out=100) - pseudo.count
    }

    if (mean(size.factors)!=1) {
        stop("size factors should be centred at unity") 
    }

    # Defining which distribution to use.
    if (dispersion==0) {
        qfun <- function(..., mean) { qpois(..., lambda=mean) }
        dfun <- function(..., mean) { dpois(..., lambda=mean) }
    } else {
        qfun <- function(..., mean) { qnbinom(..., mu=mean, size=1/dispersion) }
        dfun <- function(..., mean) { dnbinom(..., mu=mean, size=1/dispersion) }
    }
    collected.means <- collected.vars <- numeric(length(means))

    for (i in seq_along(means)) {
        m <- means[i]*size.factors
        lower <- qfun(tol, mean=m, lower.tail=TRUE)
        upper <- qfun(tol, mean=m, lower.tail=FALSE)

        # Creating a function to compute the relevant statistics.
        .getValues <- function(j) {
            ranged <- lower[j]:upper[j]
            p <- dfun(ranged, mean=m[j])
            lvals <- log2(ranged/size.factors[j] + pseudo.count)
            return(list(p=p, lvals=lvals))
        }

        # Computing the mean.
        cur.means <- numeric(length(size.factors))
        for (j in seq_along(size.factors)) { 
            out <- .getValues(j)
            cur.means[j] <- sum(out$lvals * out$p) / sum(out$p)
        }
        final.mean <- mean(cur.means)
        collected.means[i] <- final.mean
       
        # Computing the variance. Done separately to avoid
        # storing 'p' and 'lvals' in memory, but as a result
        # we need to compute these values twice.
        cur.vars <- numeric(length(size.factors))
        for (j in seq_along(size.factors)) { 
            out <- .getValues(j)
            cur.vars[j] <- sum((out$lvals - final.mean)^2 * out$p) / sum(out$p)
        }
        collected.vars[i] <- mean(cur.vars)
    }

    splinefun(collected.means, collected.vars)
}

