setGeneric("decomposeVar", function(x, fit, ...) standardGeneric("decomposeVar"))

setMethod("decomposeVar", c("matrix", "list"), function(x, fit, design=NA, subset.row=NULL)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
# last modified 6 May 2016
{
    if (is.null(design)) { 
        design <- .interceptModel(ncol(x)) 
    } else if (length(design)==1L && is.na(design)) { 
        design <- fit$design 
    }
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    QR <- .checkDesign(design)

    lout <- .Call(cxx_estimate_variance, QR$qr, QR$qraux, x, subset.row-1L)
    if (is.character(lout)) { stop(lout) }
    lmeans <- lout[[1]]
    lvar <- lout[[2]]
    
    tech.var <- fit$trend(lmeans)
    bio.var <- lvar - tech.var
    pval <- testVar(total=lvar, null=tech.var, df=nrow(design) - ncol(design))
    out <- data.frame(mean=lmeans, total=lvar, bio=bio.var, tech=tech.var,
                      p.value=pval, FDR=p.adjust(pval, method="BH"))
    rownames(out) <- rownames(x)[subset.row]
    return(out)
})

setMethod("decomposeVar", c("SCESet", "list"), function(x, fit, subset.row=NULL, ..., assay="exprs", get.spikes=FALSE) {
    out <- decomposeVar(assayDataElement(x, assay), fit, ..., subset.row=subset.row)
    if (!get.spikes) {
        nokeep <- isSpike(x, warning=FALSE)
        if (!is.null(subset.row)) { 
            nokeep <- nokeep[subset.row]
        }
        if (any(nokeep)) { 
            out$p.value[nokeep] <- NA
            out$FDR <- p.adjust(out$p.value, method="BH")
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
# last modified 6 May 2016
{
    if (missing(df)) { df <- nrow(design) - qr(design)$rank }
    if (!length(null) || !length(df)) { return(rep(NA_real_, length(total))) }
    pchisq(total/null*df, df=df, lower.tail=FALSE)
}

