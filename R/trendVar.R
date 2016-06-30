setGeneric("trendVar", function(x, ...) standardGeneric("trendVar"))

setMethod("trendVar", "matrix", function(x, trend=c("poly", "loess"), df=5, span=0.3, design=NULL, subset.row=NULL)
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 6 June 2016
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
