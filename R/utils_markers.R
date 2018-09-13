.combine_pvalues <- function(P, pval.type="any", log.p.in=FALSE, log.p.out=FALSE) 
# This function combines the p-values using Simes' method or via the IUT.
# Additional arguments are involved to specify whether the input/output should be logged.
{
    if (pval.type=="any") { 
        # Computing the Simes p-value (with NA protection).
        P <- .Call(cxx_combine_simes, P, log.p.in)
    } else {
        # Computing the IUT p-value.
        P <- P[.find_largest_col(P)]
    }

    # Deciding what to return (at this point, P is the same log-status as it was supplied).
    if (log.p.in && !log.p.out) { 
        P <- exp(P) 
    } else if (!log.p.in && log.p.out) {
        P <- log(P)
    }
    return(P)
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

