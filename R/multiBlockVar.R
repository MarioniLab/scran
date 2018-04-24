#' @export 
#' @importFrom S4Vectors DataFrame metadata
multiBlockVar <- function(x, block, trend.args=list(), dec.args=list(), ...) 
# This function fits a separate trend for each level of the blocking factor.
# The aim is to account for different mean-variance trends in a coherent manner.
# 
# written by Aaron Lun
# created 23 April 2018
{
    by.block <- split(seq_len(ncol(x)), block)
    all.out <- vector("list", length(by.block))
    names(all.out) <- names(by.block)

    for (b in names(by.block)) {
        cur.b <- by.block[[b]]
        cur.x <- x[,cur.b]

        # Estimating the technical/biological components.
        cur.fit <- do.call(trendVar, c(list(x), trend.args))
        cur.dec <- do.call(decomposeVar, c(list(x, cur.fit), dec.args))
        metadata(cur.dec)$trend <- cur.fit$trend
        
        all.out[[b]] <- cur.dec
    }

    combined <- do.call(combineVar, c(all.out, list(...)))
    stored <- do.call(DataFrame, c(lapply(all.out, I), list(check.names=FALSE)))
    combined$per.block <- stored
    return(combined)
}
