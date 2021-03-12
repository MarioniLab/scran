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

    is.empty <- groups==""
    if (any(is.empty, na.rm=TRUE)) {
        warning("replacing empty 'groups' with NA")
        groups[is.empty] <- NA
    }

    groups <- as.factor(groups)
    if (nlevels(groups) < 2L) {
        stop("need at least two unique levels in 'groups'")
    }
    groups
}

.setup_gene_names <- function(gene.names, x, subset.row) {
    if (!is.null(gene.names)) {
        .Deprecated(msg="use of custom values in 'gene.names=' is deprecated")
        if (length(gene.names)!=nrow(x)) {
            stop("length of 'gene.names' is not equal to the number of rows")
        } 
    } else {
        gene.names <- rownames(x)
        if (is.null(gene.names)) {
            gene.names <- as.character(seq_len(nrow(x)))
        }
    }

    if (!is.null(subset.row) && !is.null(gene.names)) {
        gene.names <- gene.names[.subset2index(subset.row, x, byrow=TRUE)]
    } 

    gene.names
}

#' @importFrom BiocParallel SerialParam bplapply
.pairwise_blocked_template <- function(group.vals, nblocks, direction, 
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
#' @importFrom metapod combineParallelPValues averageParallelStats
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
            com.left <- combineParallelPValues(all.left[valid.test], method="stouffer", weights=all.weight, log.p=TRUE)$p.value
            com.right <- combineParallelPValues(all.right[valid.test], method="stouffer", weights=all.weight, log.p=TRUE)$p.value

            hvt.p <- .choose_leftright_pvalues(com.left, com.right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(com.right, com.left, direction=direction)

            forward.effect <- averageParallelStats(all.forward[valid.test], all.weight)
            reverse.effect <- averageParallelStats(all.reverse[valid.test], all.weight)
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
        # The doubling can also be interpreted as a Bonferroni correction,
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

#' BH correction on log-p-values
#' 
#' Perform a Benjamini-Hochberg correction on log-transformed p-values to get log-adjusted p-values,
#' without the loss of precision from undoing and redoing the log-transformations.
#'
#' @param log.p.val Numeric vector of log-transformed p-values.
#'
#' @return A numeric vector of the same length as \code{log.p.val} containing log-transformed BH-corrected p-values.
#' @author Aaron Lun
#'
#' @examples
#' log.p.values <- log(runif(1000))
#' obs <- .logBH(log.p.values)
#' head(obs)
#'
#' ref <- log(p.adjust(exp(log.p.values), method="BH"))
#' head(ref)
#' 
#' @export
#' @rdname logBH
.logBH <- function(log.p.val) {
    o <- order(log.p.val)
    repval <- log.p.val[o] + log(length(o)/seq_along(o))
    repval <- rev(cummin(rev(repval)))
    repval[o] <- repval
    repval
}
