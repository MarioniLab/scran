testVar <- function(total, null, df, design=NULL, test=c("chisq", "f"), fit=NULL, verbose=FALSE)
# Tests that total > null given variances estimated on 'df' degrees of freedom.
# You can also give it the design matrix directly if you can't be bothered estimating 'df'.
# Obviously there's an assumption of normality here, regarding the observations from which estimation was performed.
#
# written by Aaron Lun
# created 9 February 2016
# last modified 23 December 2016
{
    if (missing(df)) { df <- nrow(design) - qr(design)$rank }
    if (!length(null) || !length(df)) { return(rep(NA_real_, length(total))) }

    test <- match.arg(test)
    if (test=="chisq") {
        if (verbose) message(sprintf("testing on %i degrees of freedom", df))
        p <- pchisq(total/null*df, df=df, lower.tail=FALSE)
    } else {
        if (is.null(fit)) { stop("fit object from trendVar() must be specified for test='f'") }
        vals <- fit$var/fit$trend(fit$mean)
        ffit <- fitFDistRobustly(vals[fit$var > 0], df1=nrow(fit$design) - ncol(fit$design))

        if (verbose) message(sprintf("testing on %i and %.2f degrees of freedom, scaled by %.2f", df, ffit$df2, ffit$scale))
        # Assumes that the scaled inverse-chisq distribution for true variances is the same.
        p <- pf(total/null/ffit$scale, df1=df, df2=ffit$df2, lower.tail=FALSE) 
    }
    return(p)
}

