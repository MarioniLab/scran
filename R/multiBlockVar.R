#' @export 
#' @importFrom S4Vectors DataFrame metadata
multiBlockVar <- function(x, block, make.tech.trend=FALSE, trend.args=list(), dec.args=list(), assay.type="logcounts", ...) 
# Fits a separate trend for each level of the blocking factor, and decomposes it.
# The aim is to account for different mean-variance trends in a coherent manner.
# 
# written by Aaron Lun
# created 23 April 2018
{
    .check_centered_SF(x, assay.type=assay.type, block)

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
        if (!make.tech.trend) {
            cur.fit <- do.call(trendVar, c(list(cur.x, assay.type=assay.type), trend.args))
        } else {
            new.trend <- do.call(makeTechTrend, c(list(x=cur.x), trend.args))
            cur.fit <- list(trend=new.trend)
        }
        cur.dec <- do.call(decomposeVar, c(list(cur.x, cur.fit, assay.type=assay.type), dec.args))
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
