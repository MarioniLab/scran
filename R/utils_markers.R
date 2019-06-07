.pairwise_blocked_template <- function(x, clusters, block=NULL, direction="any", 
    gene.names=NULL, log.p=TRUE, STATFUN, FLIPFUN, effect.name) 
{
    nblocks <- length(block)
    out.stats <- .create_output_container(clusters)

    for (i in seq_along(clusters)) {
        host <- clusters[i]
        targets <- clusters[seq_len(i-1L)]

        for (target in targets) {
            all.effect <- all.left <- all.right <- vector("list", nblocks)
            valid.test <- logical(nblocks)
            all.weight <- numeric(nblocks)

            for (b in seq_len(nblocks)) {
                out <- STATFUN(b, host, target)
                all.effect[[b]] <- out$effect
                all.weight[b] <- out$weight
                valid.test[b] <- out$valid
                all.left[[b]] <- out$left
                all.right[[b]] <- out$right
            }

            # Combining the p-values for each side across blocks.
            if (any(valid.test)) { 
                comb.args <- list(method="z", weights=all.weight[valid.test], log.p=TRUE)
                com.left <- do.call(combinePValues, c(all.left[valid.test], comb.args))
                com.right <- do.call(combinePValues, c(all.right[valid.test], comb.args))
            } else {
                com.left <- com.right <- rep(NA_real_, length(all.left[[1]]))
            }

            # Flipping left/right to get the p-value from the reversed comparison.
            hvt.p <- .choose_leftright_pvalues(com.left, com.right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(com.right, com.left, direction=direction)

            # Flipping the effect.
            com.effect <- .weighted_average_vals(all.effect, all.weight, weighted=TRUE)
            out.stats[[host]][[target]] <- .create_full_stats(effect=com.effect, p=hvt.p, 
                gene.names=gene.names, log.p=log.p, effect.name=effect.name)
            out.stats[[target]][[host]] <- .create_full_stats(effect=FLIPFUN(com.effect), p=tvh.p, 
                gene.names=gene.names, log.p=log.p, effect.name=effect.name)
        }
    }
	
    out.stats
}

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
        # The x2 can also be interpreted as a Bonferroni correction,
        # and thus is valid even if the null regions are not symmetric.
        log.p.out <- pmin(left, right) + log(2)
        return(pmin(0, log.p.out)) 
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
.create_full_stats <- function(effect, p, gene.names, log.p=TRUE, effect.name="logFC") {
    p <- as.vector(p)
    effect.list <- list(effect)
    names(effect.list) <- effect.name

    if (log.p) {
        DataFrame(effect.list, log.p.value=p, log.FDR=.logBH(p), check.names=FALSE, row.names=gene.names)
    } else {
        # Avoid underflow that might be avoided after corretion by correcting in log-space first.
        DataFrame(effect.list, p.value=exp(p), FDR=exp(.logBH(p)), check.names=FALSE, row.names=gene.names)
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
