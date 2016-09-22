setGeneric("trendVar", function(x, ...) standardGeneric("trendVar"))

setMethod("trendVar", "matrix", function(x, trend=c("loess", "semiloess"), 
    span=0.3, family="symmetric", degree=1, start=list(),
    design=NULL, subset.row=NULL)
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 11 September 2016
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    if (is.null(design)) { design <- .interceptModel(ncol(x)) }
    QR <- .checkDesign(design)

    lout <- .Call(cxx_estimate_variance, QR$qr, QR$qraux, x, subset.row - 1L)
    if (is.character(lout)) { stop(lout) }
    lmeans <- lout[[1]]
    lvar <- lout[[2]]
    names(lmeans) <- names(lvar) <- rownames(x)[subset.row]

    is.okay <- lvar > 1e-8
    kept.var <- lvar[is.okay]
    kept.means <- lmeans[is.okay]
    
    trend <- match.arg(trend)
    if (trend=="loess") { 
        kept.lvar <- log2(kept.var)
        fit <- loess(kept.lvar ~ kept.means, span=span, degree=degree, family=family)
        SUBFUN <- function(x) { 2^predict(fit, data.frame(kept.means=x)) }
    } else if (trend=="semiloess") {
        # Setting up starting parameters.
        if (is.null(start$a1)) { start$a <- 5 }
        if (is.null(start$n)) { start$n <- 5 }
        if (is.null(start$b)) { start$b <- 1 }
       
        if (length(kept.var) <= 3L) {
            stop("need at least 4 values for non-linear curve fitting")
        } 
        init.fit <- nls(kept.var ~ (a*kept.means)/(kept.means^n + b), start=start,
                        control=nls.control(warnOnly=TRUE, maxiter=500), algorithm="port",
                        lower=list(a=0, b=0, n=1))

        leftovers <- log2(kept.var/fitted(init.fit))
        loess.fit <- loess(leftovers ~ kept.means, span=span, degree=degree, family=family)
        SUBFUN <- function(x) { 
            args <- data.frame(kept.means=x)
            predict(init.fit, args) * 2^predict(loess.fit, args)
        }
    } 

    left.edge <- min(kept.means)
    right.edge <- max(kept.means)
    FUN <- function(x) {
        out <- SUBFUN(x)
        out[x < left.edge] <- SUBFUN(left.edge)
        out[x > right.edge] <- SUBFUN(right.edge)
        return(out)
    }
    return(list(mean=lmeans, var=lvar, trend=FUN, design=design))
})

.interceptModel <- function(ncells) {
    as.matrix(rep(1, ncells)) 
}

.checkDesign <- function(design) {
    if (ncol(design) >= nrow(design)) {
        stop("design matrix is not of full rank")
    }
    QR <- qr(design, LAPACK=TRUE)
    if (QR$rank!=ncol(design)) {
        stop("design matrix is not of full rank")
    }
    return(QR)
}

setMethod("trendVar", "SCESet", function(x, subset.row=NULL, ..., assay="exprs", use.spikes=TRUE) {
    if (is.null(subset.row)) {
        if (is.na(use.spikes)) {   
            subset.row <- NULL
        } else if (use.spikes) {
            subset.row <- is.spike(x)
            if (is.null(subset.row)) { 
                subset.row <- logical(nrow(x)) # no spikes at all.
            }
        } else {
            subset.row <- .spikeSubset(x, get.spikes=FALSE)
        }
    }
    out <- trendVar(assayDataElement(x, assay), ..., subset.row=subset.row)
    return(out)
})
