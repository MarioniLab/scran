#' @importFrom stats runmed
#' @export
DM <- function(mean, cv2, win.size=51) 
# Computes the distance to median for the CV2 values across all genes, 
# after fitting an abundance-dependent trend.
# 
# written by Jong Kyoung Kim
# with modifications by Aaron Lun
# created 12 March 2015 
{
    keep <- mean > 0 & !is.na(cv2) & cv2 > 0
    mean.expr <- log10(mean[keep])
    cv2.expr <- log10(cv2[keep])

    o <- order(mean.expr)
    if (win.size%%2L==0L) {
        win.size <- win.size+1L
    }
    med.trend <- runmed(cv2.expr[o], k=win.size)
    med.trend[o] <- med.trend

    dm.out <- cv2.expr - med.trend
    DM <- rep(NA_real_, length(keep))
    DM[keep] <- dm.out
    names(DM) <- names(mean)
    DM
}
 
