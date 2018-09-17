#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
pairwiseWilcox <- function(x, clusters, block=NULL, direction=c("any", "up", "down"),
    log.p=FALSE, gene.names=rownames(x), subset.row=NULL, tol=1e-8, BPPARAM=SerialParam())
# Performs pairwise Welch t-tests between clusters.
#
# written by Aaron Lun
# created 15 September 2018
{
    ncells <- ncol(x)
    clusters <- as.factor(clusters)
    if (length(clusters)!=ncells) {
        stop("length of 'clusters' does not equal 'ncol(x)'")
    }
    if (nlevels(clusters) < 2L) {
        stop("need at least two unique levels in 'clusters'")
    }

    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    if (is.null(gene.names)) {
        gene.names <- subset.row
    } else if (length(gene.names)!=nrow(x)) {
        stop("length of 'gene.names' is not equal to the number of rows")
    } else {
        gene.names <- gene.names[subset.row]
    }

    direction <- match.arg(direction)
    results <- .blocked_wilcox(x, subset.row, clusters, block=block, direction=direction, gene.names=gene.names, log.p=log.p, tol=tol, BPPARAM=BPPARAM)

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names), use.names=FALSE)
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom stats pwilcox
.blocked_wilcox <- function(x, subset.row, clusters, block=NULL, direction="any", gene.names=NULL, log.p=TRUE, tol=tol, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# t-tests between pairs of clusters within each blocking level.
{
    if (is.null(block)) {
        block <- list(`1`=seq_len(ncol(x)))
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- split(seq_along(block), block)
    }

    # Choosing the parallelization strategy.
    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row - 1L, wout)

    # Computing across blocks.
    clust.vals <- levels(clusters)
    nblocks <- length(block)
    all.stats <- all.n <- vector("list", nblocks)

    for (b in seq_along(block)) {
        chosen <- block[[b]]
        cur.clusters <- clusters[chosen]
        all.n[[b]] <- as.vector(table(cur.clusters))
        names(all.n[[b]]) <- clust.vals
        
        cur.groups <- split(chosen - 1L, cur.clusters)
        raw.stats <- bplapply(by.core, FUN=.find_overlap_exprs, x=x, by.group=cur.groups, tol=tol, BPPARAM=BPPARAM)

        cons.stats <- vector("list", length(clust.vals))
        names(cons.stats) <- clust.vals
        for (i in seq_along(cons.stats)) {
            cons.stats[[i]] <- do.call(rbind, lapply(raw.stats, "[[", i=i))
            colnames(cons.stats[[i]]) <- clust.vals[seq_len(i-1L)]
        }
        all.stats[[b]] <- cons.stats
    }

    # Organizing into some pairwise output.
    out.stats <- .create_output_container(clust.vals)
    
    for (i in seq_along(clust.vals)) {
        host <- clust.vals[i]
        targets <- clust.vals[seq_len(i-1L)]

        for (target in targets) {
            all.effect <- all.left <- all.right <- vector("list", nblocks)
            valid.test <- logical(nblocks)
            all.weight <- numeric(nblocks)

            for (b in seq_len(nblocks)) {
                cur.effect <- all.stats[[b]][[host]][,target]
                host.n <- as.double(all.n[[b]][[host]]) # numeric conversion to avoid overflow.
                target.n <- as.double(all.n[[b]][[target]])
                cur.prod <- host.n * target.n

                all.effect[[b]] <- cur.effect / cur.prod
                all.left[[b]] <- pwilcox(cur.effect, host.n, target.n, log.p=TRUE, lower.tail=TRUE)
                all.right[[b]] <- pwilcox(cur.effect, host.n, target.n, log.p=TRUE, lower.tail=FALSE)

                all.weight[b] <- cur.prod
                valid.test[b] <- host.n > 0L && target.n > 0L
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

            # Symmetric effects, hence the '1-'.
            com.effect <- .weighted_average_vals(all.effect, all.weight, weighted=TRUE)
            out.stats[[host]][[target]] <- .create_full_stats(com.effect, hvt.p, gene.names, log.p=log.p)
            out.stats[[target]][[host]] <- .create_full_stats(1-com.effect, tvh.p, gene.names, log.p=log.p)
        }
    }
	
    out.stats
}
