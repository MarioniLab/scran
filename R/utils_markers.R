.choose_leftright_pvalues <- function(left, right, direction="any")
# Choosing whether to use the p-values from the left (lower.tail), or
# the p-values from the right (upper.tail), or to make it two-sided.
# Assumes 'left' and 'right' are log-transformed p-values.
{
    if (direction=="up") {
        return(right)
    } else if (direction=="down") {
        return(left)
    } else {
        log.p.out <- pmin(left, right) + log(2)
        return(log.p.out)
    }
}

.create_output_container <- function(clust.vals) {
    out <- vector("list", length(clust.vals))
    names(out) <- clust.vals
    for (i in seq_along(clust.vals)) {
        targets <- clust.vals[-i]
        collected <- vector("list", length(targets))
        names(collected) <- targets
        host <- clust.vals[i]
        out[[host]] <- collected
    }
    return(out)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
.create_full_stats <- function(..., p, gene.names, log.p=TRUE) {
    p <- as.vector(p)
    if (log.p) {
        DataFrame(..., log.p.value=p, log.FDR=.logBH(p), check.names=FALSE, row.names=gene.names)
    } else {
        # Avoid underflow that might be avoided after corretion by correcting in log-space first.
        DataFrame(..., p.value=exp(p), FDR=exp(.logBH(p)), check.names=FALSE, row.names=gene.names)
    }
}

.logBH <- function(log.p.val) 
# Same as log(p.adjust(exp(log.p.val), method="BH")), without
# the need to undo and redo the log-transformations.
{
    o <- order(log.p.val)
    repval <- log.p.val[o] + log(length(o)/seq_along(o))
    repval <- rev(cummin(rev(repval)))
    repval[o] <- repval
    return(repval)
}

.rank_top_genes <- function(metrics) 
# This computes the rank and the minimum metric for each gene.
{
    ngenes <- nrow(metrics)
    ncon <- ncol(metrics)
    min.rank <- min.val <- rep(NA_integer_, ngenes)

    for (con in seq_len(ncon)) { 
        cur.val <- metrics[,con]
        cur.rank <- rank(cur.val, ties.method="first", na.last="keep")
        min.rank <- pmin(min.rank, cur.rank, na.rm=TRUE)
        min.val <- pmin(min.val, cur.val, na.rm=TRUE)
    }
    
    return(list(rank=min.rank, value=min.val))
}

.find_largest_col <- function(metrics) 
# Finds the indices of the largest entry in each row.
{
    metrics[is.na(metrics)] <- -Inf # protection against NAs.
    largest <- max.col(metrics)
    cbind(seq_along(largest), largest)
}

