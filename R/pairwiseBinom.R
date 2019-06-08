#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
pairwiseBinom <- function(x, clusters, block=NULL, direction=c("any", "up", "down"),
    log.p=FALSE, gene.names=rownames(x), subset.row=NULL, threshold=1e-8, BPPARAM=SerialParam())
# Performs binomial tests between clusters.
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
    results <- .blocked_binom(x, subset.row, clusters, block=block, direction=direction, 
        gene.names=gene.names, log.p=log.p, threshold=threshold, BPPARAM=BPPARAM)

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names), use.names=FALSE)
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom stats pbinom
.blocked_binom <- function(x, subset.row, clusters, block=NULL, direction="any", gene.names=NULL, log.p=TRUE, 
	threshold=1e-8, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# binomial tests between pairs of clusters within each blocking level.
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
    by.core <- .split_vector_by_workers(subset.row, wout)

    # Computing across blocks.
    clust.vals <- levels(clusters)
    nblocks <- length(block)
    all.nzero <- all.n <- vector("list", nblocks)

    for (b in seq_along(block)) {
        chosen <- block[[b]]
        cur.clusters <- clusters[chosen]
        all.n[[b]] <- as.vector(table(cur.clusters))
        names(all.n[[b]]) <- clust.vals
        
        cur.groups <- split(chosen, cur.clusters)
        raw.nzero <- bplapply(by.core, FUN=.compute_nzero_stat, x=x, by.group=cur.groups, threshold=threshold, BPPARAM=BPPARAM)
        cons.nzero <- do.call(rbind, raw.nzero)
        colnames(cons.nzero) <- clust.vals
        all.nzero[[b]] <- cons.nzero
    }

    # Not exactly equal to binom.test(); the two-sided p-value from binom.test()
    # cannot be performed by any combination of the one-sided p-values. This 
    # makes it impossible to behave with directional Stouffer's method in 
    # .pairwise_blocking_template(), and just generally gums up the works.
    STATFUN <- function(b, host, target) {    
        host.nzero <- all.nzero[[b]][,host]
        target.nzero <- all.nzero[[b]][,target]
        host.n <- all.n[[b]][[host]]
        target.n <- all.n[[b]][[target]]

        size <- host.nzero + target.nzero
        p <- host.n/(host.n + target.n)

        list(
            effect=unname(log2((host.nzero + 0.5)/(host.n + 0.5) * (target.n + 0.5)/(target.nzero + 0.5))),
            weight=as.double(host.n)*as.double(target.n),
            valid=host.n > 0L && target.n > 0L,
            left=pbinom(host.nzero, size, p, log.p=TRUE),
            right=pbinom(host.nzero - 1, size, p, lower.tail=FALSE, log.p=TRUE)
        )
    }

    .pairwise_blocked_template(x, clust.vals, nblocks=length(block), direction=direction, 
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, FLIPFUN=function(x) -x, effect="logOR")
}

#' @importFrom scater nexprs
.compute_nzero_stat <- function(x, by.group, rows, threshold) {
    collected <- lapply(by.group, FUN=function(s) {
        nexprs(x, subset_col=s, subset_row=rows, byrow=TRUE, detection_limit=threshold)
    })
    do.call(cbind, collected)
}
