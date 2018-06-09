#' @export 
#' @importFrom S4Vectors DataFrame metadata
multiBlockVar <- function(x, block, trend.args=list(), dec.args=list(), assay.type="logcounts", ...) 
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
        cur.fit <- do.call(trendVar, c(list(cur.x, assay.type=assay.type), trend.args))
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

#' @export
#' @importFrom scater librarySizeFactors
#' @importFrom SingleCellExperiment sizeFactorNames
#' @importFrom BiocGenerics sizeFactors sizeFactors<- normalize
#' @importFrom stats ave
multiBlockNorm <- function(x, block, ...) 
# Adjusts size factors so that the spike-in size factors have
# the same mean as the reference size factors _within each batch_.
# Ensures that the abundances are comparable in multiBlockVar.
#
# written by Aaron Lun
# created 4 June 2018
{
    ref <- sizeFactors(x)
    if (is.null(ref)) {
        ref <- librarySizeFactors(x)
        warning("using library sizes as size factors")
    } else {
        # Centering just in case.
        ref <- ref / mean(ref)
        sizeFactors(x) <- ref
    }
    ref.mean <- ave(ref, block, FUN=mean)

    for (sf in sizeFactorNames(x)) {
        current <- sizeFactors(x, sf)
        cur.mean <- ave(current, block, FUN=mean)
        sizeFactors(x, sf) <- current * ref.mean / cur.mean
    }
   
    x <- normalize(x, ..., centre_size_factors=FALSE)
    return(x)
}
