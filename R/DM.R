quantile.trend <- function(y, x = NULL, nsplits=NULL, size=20, overlap=NULL, distance=NULL, quant=0.95,
                          alignment="center", FUN=function(x) { quantile(x, quant) }, ...)
# Fits a quantile loess curve to a response 'y' with respect to a covariate 'x'.
#
# Taken from source("http://www.r-statistics.com/wp-content/uploads/2010/04/Quantile.loess_.r.txt")
# with modifications by Aaron Lun
# last modified 13 January 2016
{
    if(!is.null(nsplits)) {
        size <- ceiling(length(y)/nsplits)
    }
    if (is.null(distance)) { 
        distance <- size
    }
    if (!is.null(overlap)) { 
        distance <- size * (1-overlap)
    }
    if(is.null(x)) { 
        x <- index(y)
    } 
    zoo.y <- zoo(x = y, order.by = x)
    new.y <- rollapply(zoo.y, width = size, FUN=FUN, by=distance, align=alignment) 
    new.x <- attributes(new.y)$index  
    return(list(y = new.y, x = new.x))
}

DM <- function(mean, cv2, win.size=50) 
# Computes the distance to median for the CV2 values across all genes, 
# after fitting an abundance-dependent trend.
# 
# written by Jong Kyoung Kim
# with modifications by Aaron Lun
# created 12 March 2015 
# last modified 13 January 2016
{
    keep <- mean > 0 & !is.na(cv2) & cv2 > 0
    mean.expr <- log10(mean[keep])
    cv2.expr <- log10(cv2[keep])
    qloess <- quantile.trend(cv2.expr, mean.expr, quant=.5, size=win.size, overlap=.5)
    dm.out <- cv2.expr - approx(x=qloess$x, y=qloess$y, xout=mean.expr)$y
    DM <- rep(NA_real_, length(keep))
    DM[keep] <- dm.out
    names(DM) <- names(mean)
    DM
}
 
