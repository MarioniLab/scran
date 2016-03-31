setGeneric("decomposeVar", function(x, fit, ...) standardGeneric("decomposeVar"))

setMethod("decomposeVar", c("ANY", "list"), function(x, fit, design=NA)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
# last modified 21 February 2016
{
    x <- as.matrix(x)
    if (is.null(design)) { design <- .interceptModel(ncol(x)) }
    else if (length(design)==1L && is.na(design)) { design <- fit$design }
    lmeans <- rowMeans(x)
    lvar <- .estimateVariance(design, x)
    tech.var <- fit$trend(lmeans)
    bio.var <- lvar - tech.var
    return(data.frame(mean=lmeans, total=lvar, bio=bio.var, tech=tech.var))
})

setMethod("decomposeVar", c("SCESet", "list"), function(x, fit, ..., assay="exprs", get.spikes=FALSE) {
    cur.assay <- .getUsedMatrix(x, assay, get.spikes=TRUE)
    out <- decomposeVar(cur.assay, fit, ...)
    if (!get.spikes) {
        nokeep <- is.spike(x)
        if (any(nokeep)) { 
            out[nokeep,] <- NA
        }
    }
    return(out)
})

testVar <- function(total, null, df, design=NULL) 
# Tests that total > null given variances estimated on 'df' degrees of freedom.
# You can also give it the design matrix directly if you can't be bothered estimating 'df'.
# Obviously there's an assumption of normality here, regarding the observations from which estimation was performed.
#
# written by Aaron Lun
# created 9 February 2016    
{
    if (missing(df)) { df <- nrow(design) - qr(design)$rank }
    if (!length(null) || !length(df)) { return(rep(NA_real_, length(total))) }
    pchisq(total/null*df, df=df, lower.tail=FALSE)
}

