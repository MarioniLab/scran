#' @importFrom stats nls loess predict fitted
#' @importFrom limma fitFDistRobustly
#' @importFrom BiocParallel SerialParam
.trend_var <- function(x, method=c("loess", "spline"), parametric=FALSE, 
                       loess.args=list(), spline.args=list(), nls.args=list(),
                       block=NULL, design=NULL, weighted=TRUE, 
                       min.mean=0.1, subset.row=NULL, BPPARAM=SerialParam())
# Fits a smooth trend to the technical variability of the log-CPMs
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
{
    stats.out <- .get_var_stats(x, block=block, design=design, subset.row=subset.row, BPPARAM=BPPARAM)

    # Filtering out zero-variance and low-abundance genes.
    is.okay <- !is.na(stats.out$vars) & stats.out$vars > 1e-8 & stats.out$means >= min.mean 
    kept.vars <- stats.out$vars[is.okay]
    kept.means <- stats.out$means[is.okay]

    # Figuring out the d.f. to keep for weighting.
    if (is.matrix(is.okay)) { 
        kept.resid <- stats.out$resid.df[which(is.okay, arr.ind=TRUE)[,2]]
    } else {
        kept.resid <- rep(stats.out$resid.df, sum(is.okay))
    }

    # Fitting a parametric curve to try to flatten the shape.
    # This is of the form y = a*x/(x^n + b), but each coefficent is actually set
    # to exp(*) to avoid needing to set lower bounds.
    if (parametric) { 
        if (length(kept.vars) <= 3L) {
            stop("need at least 4 values for non-linear curve fitting")
        } 

        nls.args <- .setup_nls_args(nls.args, start.args=list(vars=kept.vars, means=kept.means))
        nls.args$formula <- kept.vars ~ (exp(A)*kept.means)/(kept.means^(1+exp(N)) + exp(B))
        if (weighted) {
            nls.args$weights <- kept.resid
        }

        init.fit <- do.call(nls, nls.args)
        to.fit <- log(kept.vars/fitted(init.fit))
        SUBSUBFUN <- function(x) { predict(init.fit, data.frame(kept.means=x)) }
    } else {
        if (length(kept.vars) < 2L) {
            stop("need at least 2 values for non-parametric curve fitting")
        } 
        to.fit <- log(kept.vars)
        left.edge <- min(kept.means)
        SUBSUBFUN <- function(x) { pmin(1, x/left.edge) } # To get a gradient from 0 to 1 below the supported range.
    } 
    
    # Fitting loess or splines to the remainder. Note that we need kept.means as a variable in the formula for predict() to work!
    method <- match.arg(method)
    if (method=="loess") { 
        loess.args <- .setup_loess_args(loess.args)
        loess.args$formula <- to.fit ~ kept.means 
        if (weighted) {
            loess.args$weights <- kept.resid
        }

        after.fit <- do.call(loess, loess.args)
        PREDICTOR <- function(x) { predict(after.fit, data.frame(kept.means=x)) }
    } else {
        spline.args <- .setup_spline_args(spline.args)
        spline.args$x <- kept.means
        spline.args$y <- to.fit
        if (weighted) {
            spline.args$w <- kept.resid/max(kept.resid) # For some reason, robustSmoothSpline only accepts weights in [0, 1].
        }

        after.fit <- do.call(aroma.light::robustSmoothSpline, spline.args)
        PREDICTOR <- function(x) { predict(after.fit, x)$y }
    }

    # Only trusting the parametric SUBSUBFUN for extrapolation; restricting non-parametric forms within the supported range.
    left.edge <- min(kept.means)
    right.edge <- max(kept.means)
    SUBFUN <- function(x) { 
        both.bounded <- pmax(pmin(x, right.edge), left.edge)
        exp(PREDICTOR(both.bounded)) * SUBSUBFUN(x)
    }

    # Estimating the df2, as well as scale shift from estimating mean of logs (assuming shape of trend is correct).
    leftovers <- kept.vars/SUBFUN(kept.means)
    f.fit <- fitFDistRobustly(leftovers, df1=kept.resid)
    f.df2 <- f.fit$df2
    
    # We don't just want the scaling factor, we want the scaled mean of the F-distribution (see explanation below).
    if (is.infinite(f.df2)) { 
        f.scale <- f.fit$scale
    } else if (f.df2 > 2) {
        f.scale <- f.fit$scale * f.df2/(f.df2 - 2) 
    } else {
        warning("undefined expectation for small df2")
        f.scale <- f.fit$scale
    }

    FUN <- function(x) { 
        output <- SUBFUN(x) * f.scale
        names(output) <- names(x)
        return(output)
    }
    stats.out$trend <- FUN
    stats.out$df2 <- f.df2
    return(stats.out)
}

#########################################################
# Computing variance and mean based on blocking factors #
#########################################################

#' @importFrom BiocParallel bplapply
.get_var_stats <- function(x, block, design, subset.row, BPPARAM) { 
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- wout
    for (core in seq_along(wout)) {
        by.core[[core]] <- subset.row[wout[[core]]]
    }
    recorder <- list()

    if (!is.null(block)) { 
        if (ncol(x)!=length(block)) {
            stop("length of 'block' should be the same as 'ncol(x)'")
        }
        recorder$block <- block

        # Checking residual d.f.
        by.block <- split(seq_len(ncol(x))-1L, block, drop=TRUE)
        resid.df <- lengths(by.block) - 1L
        if (all(resid.df<=0L)){ 
            stop("no residual d.f. in any level of 'block' for variance estimation")
        }

        # Calculating the statistics for each block. 
        raw.stats <- bplapply(by.core, FUN=.fit_oneway, by.block=by.block, x=x, BPPARAM=BPPARAM)
        means <- do.call(rbind, lapply(raw.stats, FUN="[[", i=1))
        vars <- do.call(rbind, lapply(raw.stats, FUN="[[", i=2))
        dimnames(means) <- dimnames(vars) <- list(rownames(x)[subset.row], names(by.block))

    } else if (!is.null(design)) {
        checked <- .make_var_defaults(x, fit=NULL, design=design)
        design <- checked$design
        recorder$design <- design

        # Checking residual d.f.
        resid.df <- nrow(design) - ncol(design)
        if (resid.df <= 0L) {
            stop("no residual d.f. in 'design' for variance estimation") 
        }
       
        # Calculating the residual variance of the fitted linear model. 
        QR <- .ranksafe_qr(design)
       
        raw.stats <- bplapply(by.core, FUN=.fit_linear_model, qr=QR$qr, qraux=QR$qraux, x=x, BPPARAM=BPPARAM)
        means <- unlist(lapply(raw.stats, FUN="[[", i=1))
        vars <- unlist(lapply(raw.stats, FUN="[[", i=2))
        names(means) <- names(vars) <- rownames(x)[subset.row]

    } else {
        resid.df <- ncol(x) - 1L
        if (resid.df <= 0L) {
            stop("no residual d.f. in 'x' for variance estimation")
        }

        by.block <- list(seq_len(ncol(x))-1L)
        raw.stats <- bplapply(by.core, FUN=.fit_oneway, by.block=by.block, x=x, BPPARAM=BPPARAM)
        means <- unlist(lapply(raw.stats, FUN="[[", i=1))
        vars <- unlist(lapply(raw.stats, FUN="[[", i=2))
        names(means) <- names(vars) <- rownames(x)[subset.row]

    }

    return(c(list(vars=vars, means=means, resid.df=resid.df), recorder))
}

# Helper functions to ensure that the scran namespace is carried into bplapply.
.fit_oneway <- function(by.block, x, chosen) {
    .Call(cxx_fit_oneway, by.block, x, chosen - 1L)
}

.fit_linear_model <- function(qr, qraux, x, chosen) {
    .Call(cxx_fit_linear_model, qr, qraux, x, chosen - 1L, FALSE)
}

#########################################################
# Computing NLS starting points for parametric fitting. #
#########################################################

#' @importFrom stats coef lm fitted
.get_nls_starts <- function(vars, means, grad.prop=0.5, grid.length=100, grid.max=10) {
    lvars <- log2(vars)
    fit <- loess(lvars ~ means, degree=1)

    # Getting a rough peak location from the fitted values.
    toppt <- which.max(fitted(fit))
    top.x <- means[toppt]
    top.y <- 2^fitted(fit)[toppt]

    # Getting the initial gradient, from all points halfway to the peak.
    min.x <- min(means)
    lower.set <- means < min.x + (top.x - min.x)*grad.prop
    lower.vars <- vars[lower.set]
    lower.means <- means[lower.set]
    line.fit <- lm(lower.vars ~ lower.means)
    grad <- coef(line.fit)[2]

    # If 'n' is known, the others can be solved by knowing the gradient near zero (grad=a/b), 
    # knowing the location of the peak (top.x^n=b/(n-1)) and the height of the peak (top.y=a*top.x*(n-1)/nb).
    # The last could be used to solve for 'n', but this seems unstable, so we do a grid search instead.
    getB <- function(n) (n-1)*top.x^n
    getA <- function(n) getB(n) * grad

    # Getting a rough 'n' estimate with a grid search.
    grid.pts <- seq(from=0, to=grid.max, length.out=grid.length)
    grid.ss <- vapply(grid.pts, FUN=function(i) {
        n <- 2^i
        resid <- vars - (getA(n)*means)/(means^n + getB(n))
        sum(resid^2)
    }, FUN.VALUE=0)

    n <- 2^grid.pts[which.min(grid.ss)]
    b <- unname(getB(n))
    a <- unname(getA(n))
    list(n=n, b=b, a=a)
}

##############################
# Setting default arguments. #
##############################

#' @importFrom stats loess
.setup_loess_args <- function(loess.args) {
    if (is.null(loess.args)) { 
        return(list())
    }

    loess.call <- do.call(call, c("loess", loess.args))
    loess.call <- match.call(loess, loess.call)
    loess.args <- as.list(loess.call)[-1] # expand to full argument names.

    altogether <- c(loess.args, list(span=0.3, family="symmetric", degree=1))
    keep <- !duplicated(names(altogether)) # loess.args are favoured.
    return(altogether[keep])
}

.setup_spline_args <- function(spline.args) {
    if (is.null(spline.args)) { 
        return(list())
    }

    spline.call <- do.call(call, c("robustSmoothSpline", spline.args))
    spline.call <- match.call(aroma.light::robustSmoothSpline, spline.call)
    spline.args <- as.list(spline.call)[-1]

    altogether <- c(spline.args, list(df=4, method="symmetric"))
    keep <- !duplicated(names(altogether))  # spline.args are favoured.
    return(altogether[keep])
}

#' @importFrom stats nls nls.control
.setup_nls_args <- function(nls.args, start.args) {
    if (is.null(nls.args)) { 
        return(list())
    }

    nls.call <- do.call(call, c("nls", nls.args))
    nls.call <- match.call(nls, nls.call)
    nls.args <- as.list(nls.call)[-1]

    control <- nls.control(warnOnly=TRUE, maxiter=500)
    raw.start <- do.call(.get_nls_starts, start.args)
    start <- list(A=log(raw.start$a), 
                  B=log(raw.start$b), 
                  N=log(pmax(1e-8, raw.start$n-1))) # reflects enforced positivity in formula.

    altogether <- c(nls.args, list(control=control, start=start))
    keep <- !duplicated(names(altogether)) # nls.args are favoured.
    return(altogether[keep])
}

#########################
# Setting up S4 methods #
#########################

#' @export
setGeneric("trendVar", function(x, ...) standardGeneric("trendVar"))

#' @export
setMethod("trendVar", "ANY", .trend_var)

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment isSpike
setMethod("trendVar", "SingleCellExperiment", function(x, subset.row=NULL, ..., assay.type="logcounts", use.spikes=TRUE) {
    .check_centered_SF(x, assay.type=assay.type)
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Can use only spikes, everything but spikes, or everything. 
    if (!is.na(use.spikes)) {
        is.spike <- isSpike(x)
        if (is.null(is.spike)) {
            is.spike <- logical(nrow(x))
        }
        is.spike <- which(is.spike)
        if (use.spikes) {
            if (!length(is.spike)) {
                stop("no spike-in transcripts present for 'use.spikes=TRUE'")
            }
            subset.row <- intersect(subset.row, is.spike)                
        } else {
            subset.row <- setdiff(subset.row, is.spike)
        }
    }

    out <- .trend_var(assay(x, i=assay.type), ..., subset.row=subset.row)
    return(out)
})

# EXPLANATION OF TREND:
# To wit, the technical variance is modelled as:
#    var ~ s0*F(df1, df2)
# When we log it, we get:
#    log(var) ~ log(s0) + log(F(df1, df2))
# When we fit a trend to the log-variances, we get E(log(var)) for any mu.
# This can also be expressed as:
#    E(log(var)) = log(s0) + E(log(F(df1, df2)))
# So, dividing var by exp[E(log(var))] would give us:
#    var/exp[E(log(var))] ~ F(df1, df2)/exp[E(log(F(df1, df2)))]
# The estimated scale from fitFDistRobustly represents 1/exp[E(log(F(df1, df2)))].
# Thus, if we want E(var), we should be getting:
#    E(var) = s0 * df2/(df2 - 2)
#           = exp[E(log(var))] * 1/exp[E(log(F(df1, df2)))] * df2/(df2 - 2)
# ... the first term of which is 2^predict, and the product of the latter two terms is 'f.scale'.

