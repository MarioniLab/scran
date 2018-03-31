#' @importFrom stats pchisq pf
#' @export
testVar <- function(total, null, df, design=NULL, test=c("chisq", "f"), second.df=NULL, log.p=FALSE)
# Tests that total > null given variances estimated on 'df' degrees of freedom.
# You can also give it the design matrix directly if you can't be bothered estimating 'df'.
# Obviously there's an assumption of normality here, regarding the observations from which estimation was performed.
#
# written by Aaron Lun
# created 9 February 2016
{
    if (missing(df)) { 
        df <- nrow(design) - qr(design)$rank 
    }

    test <- match.arg(test)
    if (test=="chisq") {
        p <- pchisq(total/null*df, df=df, lower.tail=FALSE, log.p=log.p)
        p[total==0 & null==0] <- 1
        p[rep(df==0, length.out=length(p))] <- NA_real_
    } else {
        if (is.null(second.df)) { 
            stop("second df from trendVar() must be specified for test='f'") 
        }

        # remove mean effect, so just left with scaling factor for F-distribution. 
        if (!is.infinite(second.df)) { 
            mean.adj <- second.df/(second.df-2)
        } else {
            mean.adj <- 1
        }
        scaling <- null/mean.adj 

        # Assumes that the scaled inverse-chisq distribution for true variances is the same for incoming genes.
        p <- pf(total/scaling, df1=df, df2=second.df, lower.tail=FALSE, log.p=log.p) 
        p[total==0 & null==0] <- 1
        p[rep(df==0 | second.df==0, length.out=length(p))] <- NA_real_
    }

    return(p)
}

