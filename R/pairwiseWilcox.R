#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
pairwiseWilcox <- function(x, clusters, block=NULL, direction=c("any", "up", "down"),
    log.p=FALSE, gene.names=rownames(x), subset.row=NULL, tol=1e-8, BPPARAM=SerialParam())
# Performs pairwise Wilcoxon rank sum tests between clusters.
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
    results <- .blocked_wilcox(x, subset.row, clusters, block=block, direction=direction, 
        gene.names=gene.names, log.p=log.p, tol=tol, BPPARAM=BPPARAM)

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names), use.names=FALSE)
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom stats pnorm 
.blocked_wilcox <- function(x, subset.row, clusters, block=NULL, direction="any", gene.names=NULL, 
    log.p=TRUE, tol=1e-8, BPPARAM=SerialParam())
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
    all.stats <- all.ties <- all.n <- vector("list", nblocks)

    for (b in seq_along(block)) {
        chosen <- block[[b]]
        cur.clusters <- clusters[chosen]
        all.n[[b]] <- as.vector(table(cur.clusters))
        names(all.n[[b]]) <- clust.vals
        
        cur.groups <- split(chosen - 1L, cur.clusters)
        bpl.out <- bplapply(by.core, FUN=.find_overlap_exprs, x=x, by.group=cur.groups, tol=tol, BPPARAM=BPPARAM)
        raw.stats <- lapply(bpl.out, "[[", i=1)
        raw.ties <- lapply(bpl.out, "[[", i=2)

        cons.stats <- cons.ties <- vector("list", length(clust.vals))
        names(cons.stats) <- names(cons.ties) <- clust.vals
        for (i in seq_along(cons.stats)) {
            cons.stats[[i]] <- do.call(rbind, lapply(raw.stats, "[[", i=i))
            cons.ties[[i]] <- do.call(rbind, lapply(raw.ties, "[[", i=i))
            colnames(cons.stats[[i]]) <- colnames(cons.ties[[i]]) <- clust.vals[seq_len(i-1L)]
        }

        all.stats[[b]] <- cons.stats
        all.ties[[b]] <- cons.ties
    }

    # This looks at every level of the blocking factor and performs
    # Wilcoxon tests between pairs of clusters within each blocking level.
    STATFUN <- function(b, host, target) {
        cur.effect <- all.stats[[b]][[host]][,target]
        host.n <- as.double(all.n[[b]][[host]]) # numeric conversion to avoid overflow.
        target.n <- as.double(all.n[[b]][[target]])
        cur.prod <- host.n * target.n

        output <- list(effect=cur.effect / cur.prod, weight=cur.prod)

        # 'cur.prod' is still nominally integer; we use 1.5 to avoid numerical imprecision upon an exact comparison.
        output$valid <- cur.prod > 1.5

        # Approximate Wilcoxon with continuity correction: ripped straight from wilcox.test() in stats.
        TIESUM <- all.ties[[b]][[host]][,target]
        z <- cur.effect - cur.prod/2
        SIGMA <- sqrt((host.n * target.n/12) * ((host.n + target.n + 1) - TIESUM/((host.n + target.n) * (host.n + target.n - 1))))

        # using 0.25 to avoid numerical imprecision; z should go up in units of 0.5's.
        CORRECTION <- if (direction=="any") ifelse(abs(z) < 0.25, 0, 0.5) else 0.5 
        output$left <- pnorm((z + CORRECTION)/SIGMA, log.p=TRUE)
        output$right <- pnorm((z - CORRECTION)/SIGMA, log.p=TRUE, lower.tail=FALSE)

        output
    }

    .pairwise_blocked_template(x, clust.vals, nblocks=length(block), direction=direction, 
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, FLIPFUN=function(x) 1-x, effect.name="overlap")
}

.find_overlap_exprs <- function(x, subset.row, by.group, tol) 
# Pass all arguments explicitly rather than via function environment
# (preserve scran namespace, avoid duplication of memory in bplapply).
{
    .Call(cxx_overlap_exprs, x, subset.row, by.group, tol)
}
