.setup_groups <- function(groups, x, restrict, exclude) { 
    ncells <- ncol(x)
    if (length(groups)!=ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }

    # Only dropping levels if 'restrict' or 'exclude' are specified,
    # to respect any empty levels in the input grouping.
    if (!is.null(restrict)) {
        groups[!groups%in% restrict] <- NA
        if (is.factor(groups)) groups <- droplevels(groups)
    }
    if (!is.null(exclude)) {
        groups[groups %in% exclude] <- NA
        if (is.factor(groups)) groups <- droplevels(groups)
    }

    groups <- as.factor(groups)
    if (nlevels(groups) < 2L) {
        stop("need at least two unique levels in 'groups'")
    }
    groups
}

.setup_gene_names <- function(gene.names, x, subset.row) {
    if (is.null(gene.names)) {
        subset.row
    } else if (length(gene.names)!=nrow(x)) {
        stop("length of 'gene.names' is not equal to the number of rows")
    } else if (!is.null(subset.row)) {
        gene.names[subset.row]
    } else {
        gene.names
    }
}

#' @importFrom BiocParallel SerialParam bplapply
.pairwise_blocked_template <- function(x, group.vals, nblocks, direction, 
    gene.names, log.p, STATFUN, effect.name, BPPARAM=SerialParam()) 
{
    bp.out <- bplapply(seq_along(group.vals), FUN=.pairwise_blocked_internal,
        group.vals=group.vals, nblocks=nblocks, direction=direction,
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, 
        effect.name=effect.name, BPPARAM=BPPARAM)

    output <- list(
        statistics=unlist(lapply(bp.out, "[[", i=1), recursive=FALSE),
        pairs=do.call(rbind, lapply(bp.out, "[[", i=2))
    )
    .reorder_pairwise_output(output, group.vals)
}

#' @importFrom S4Vectors DataFrame
.pairwise_blocked_internal <- function(i, group.vals, nblocks, direction="any", 
    gene.names=NULL, log.p=TRUE, STATFUN, effect.name) 
{
    host <- group.vals[i]
    targets <- group.vals[seq_len(i-1L)]
    collected.stats <- collected.pairs <- list()
    counter <- 1L

    for (target in targets) {
        all.forward <- all.reverse <- all.left <- all.right <- vector("list", nblocks)
        valid.test <- logical(nblocks)
        all.weight <- numeric(nblocks)

        for (b in seq_len(nblocks)) {
            out <- STATFUN(b, host, target)
            all.forward[[b]] <- out$forward
            all.reverse[[b]] <- out$reverse
            all.weight[b] <- out$weight
            valid.test[b] <- out$valid
            all.left[[b]] <- out$left
            all.right[[b]] <- out$right
        }

        # Combining the p-values for each side across blocks.
        if (any(valid.test)) { 
            all.weight <- all.weight[valid.test]
            comb.args <- list(method="z", weights=all.weight, log.p=TRUE)
            com.left <- do.call(combinePValues, c(all.left[valid.test], comb.args))
            com.right <- do.call(combinePValues, c(all.right[valid.test], comb.args))

            hvt.p <- .choose_leftright_pvalues(com.left, com.right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(com.right, com.left, direction=direction)

            forward.effect <- .weighted_average_vals(all.forward[valid.test], all.weight)
            reverse.effect <- .weighted_average_vals(all.reverse[valid.test], all.weight)
        } else {
            hvt.p <- tvh.p <- forward.effect <- reverse.effect <- rep(NA_real_, length(all.left[[1]]))
            warning(paste("no within-block comparison between", host, "and", target))
        }

        collected.stats[[counter]] <- list(
            .create_full_stats(effect=forward.effect, p=hvt.p, gene.names=gene.names, log.p=log.p, effect.name=effect.name),
            .create_full_stats(effect=reverse.effect, p=tvh.p, gene.names=gene.names, log.p=log.p, effect.name=effect.name)
        )
        collected.pairs[[counter]] <- DataFrame(first=c(host, target), second=c(target, host))
        counter <- counter + 1L
    }

    list(
        statistics=unlist(collected.stats, recursive=FALSE), 
        pairs=do.call(rbind, collected.pairs)
    )
}

.reorder_pairwise_output <- function(output, levels) {
    o <- order(
        match(output$pairs$first, levels), 
        match(output$pairs$second, levels)
    )
    output$statistics <- output$statistics[o]
    output$pairs <- output$pairs[o,]
    output
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

.weighted_average_vals <- function(vals, weights) {
    combined <- total.weight <- 0
    for (x in seq_along(vals)) {
        cur.weights <- weights[[x]]
        product <- vals[[x]] * cur.weights

        # avoid problems with NA values that have zero weight.
        product[is.na(product) & cur.weights==0] <- 0

        combined <- combined + product
        total.weight <- total.weight + cur.weights
    }
    combined/total.weight
}
