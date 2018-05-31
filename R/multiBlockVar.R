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
        if (length(cur.b) < 2L) {
            warning("fewer than two cells in a block")
            next
        }
        cur.x <- x[,cur.b]

        # Estimating the technical/biological components.
        cur.fit <- do.call(trendVar, c(list(cur.x), trend.args))
        cur.dec <- do.call(decomposeVar, c(list(cur.x, cur.fit), dec.args))
        metadata(cur.dec)$trend <- cur.fit$trend
        
        all.out[[b]] <- cur.dec
    }
    
    # Getting rid of empty blocks.
    is.empty <- vapply(all.out, FUN=is.null, FUN.VALUE=TRUE)
    if (!any(!is.empty)) {
        stop("no block contains enough cells for variance estimation")
    }
    all.out <- all.out[!is.empty]

    combined <- do.call(combineVar, c(all.out, list(...)))
    stored <- do.call(DataFrame, c(lapply(all.out, I), list(check.names=FALSE)))
    combined$per.block <- stored
    return(combined)
}
