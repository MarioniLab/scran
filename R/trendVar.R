setGeneric("trendVar", function(x, ...) standardGeneric("trendVar"))

setMethod("trendVar", "matrix", function(x, trend=c("poly", "loess"), df=5, span=0.3, design=NULL)
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 21 February 2016
{
    lmeans <- rowMeans(x)
    if (is.null(design)) { design <- .interceptModel(ncol(x)) }
    lvar <- .estimateVariance(design, x)

    is.okay <- lvar > 1e-8
    kept.means <- lmeans[is.okay]
    llvar <- log2(lvar)[is.okay]
    
    trend <- match.arg(trend)
    if (trend=="loess") { 
        fit <- loess(llvar ~ kept.means, span=span, degree=1, family="symmetric")
    } else if (trend=="poly") {
        fit <- lm(llvar ~ poly(kept.means, df=df))
    } 

    left.edge <- which.min(kept.means)
    right.edge <- which.max(kept.means)
    FUN <- function(x) {
        out <- predict(fit, data.frame(kept.means=x))
        out[x < kept.means[left.edge]] <- fitted(fit)[left.edge]
        out[x > kept.means[right.edge]] <- fitted(fit)[right.edge]
        return(2^out)
    }
    return(list(mean=lmeans, var=lvar, trend=FUN, design=design))
})

.interceptModel <- function(ncells) {
    as.matrix(rep(1, ncells)) 
}

setMethod("trendVar", "SCESet", function(x, ..., assay="exprs", use.spikes=TRUE) {
    if (is.na(use.spikes)) {   
        cur.assay <- assayDataElement(x, assay)
    } else if (use.spikes) {
        cur.assay <- spikes(x, assay)
    } else {
        cur.assay <- .getUsedMatrix(x, assay)
    }
    out <- trendVar(cur.assay, ...)
    return(out)
})

.estimateVariance <- function(X, y) {
    if (!nrow(y)) { return(numeric(0)) }
    else if (!ncol(y)) { return(rep(NaN, nrow(y))) }
    fit <- lm.fit(x=X, y=t(y))
    return(colMeans(fit$effects[-fit$qr$pivot[seq_len(fit$rank)],,drop=FALSE]^2))
}

