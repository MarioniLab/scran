#' Score marker genes
#'
#' Compute various summary scores for potential marker genes to distinguish between groups of cells.
#'
#' @param x A matrix-like object containing log-normalized expression values, with genes in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix in its assays.
#' @param groups A factor or vector containing the identity of the group for each cell in \code{x}.
#' @param block,design Further arguments to be passed to \code{\link{pairwiseTTest}} and related functions.
#' @param weight.fun Function indicating how the statistics from different comparisons should be weighted.
#' This should accept a vector of integers containing the number of cells involved in the comparison,
#' and return a vector of equal length containing the weights. 
#' @param row.data A DataFrame of length equal to \code{nrow(x)}, containing extra information to insert into each DataFrame.
#' @param full.stats Logical scalar indicating whether the statistics from the pairwise comparisons should be directly returned.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the calculations should be parallelized.
#' @param ... For the generic, further arguments to pass to individual methods.
#' 
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#'
#' @return
#' A List of DataFrames containing marker scores for each gene in each group.
#' Each DataFrame corresponds to a group and each row corresponds to a gene in \code{x}.
#' See Details for information about the individual columns.
#'
#' @details
#' Compared to \code{\link{findMarkers}}, this function represents a simpler and more intuitive summary of the differences between the groups.
#' We do this by realizing that the p-values for these types of comparisons are largely meaningless;
#' individual cells are not meaningful units of experimental replication, while the groups themselves are defined from the data.
#' Thus, by discarding the p-values, we can simplify our marker selection by focusing only on the effect sizes between groups.
#' 
#' Here, the strategy is to perform pairwise comparisons between each pair of groups to obtain various effect sizes.
#' For each group \eqn{X}, we summarize the effect sizes across all pairwise comparisons involving that group, e.g., mean, min, max and so on.
#' This yields a DataFrame for each group where each column contains a different summarized effect and each row corresponds to a gene in \code{x}.
#' Reordering the rows by the summary of choice can yield a ranking of potential marker genes for downstream analyses.
#' 
#' In addition, we report the mean log-expression of all cells in \eqn{X}, as well as the grand mean of mean log-expression values for all other groups.
#' This can be used to easily compute an overall log-fold change though ranking is best performed on one of the effect sizes described below.
#' We also report the proportion of cells with detectable expression in \eqn{X} and the mean proportion for all other groups.
#'
#' @section Choice of effect sizes:
#' The \code{logFC.cohen} columns contain the standardized log-fold change, i.e., Cohen's d.
#' For each pairwise comparison, this is defined as the difference in the mean log-expression for each group scaled by the root of the pooled variance across the groups.
#' The standardization is analogous to the calculation of the t-statistic and avoids spuriously large effect sizes from highly variable groups.
#' We can also interpret Cohen's d as the number of standard deviations between the two group means.
#' 
#' The \code{AUC} columns contain the area under the curve.
#' This is the probability that a randomly chosen observation in one group is greater than a randomly chosen observation in the other group.
#' The AUC is closely related to the U-statistic used in the Wilcoxon rank sum test.
#' Values greater than 0.5 indicate that a gene is upregulated in the first group.
#'
#' The key difference between the AUC and Cohen's d is that the former is less sensitive to the variance within each group.
#' The clearest example is that of two distributions that exhibit no overlap, where the AUC is the same regardless of the variance of each distribution.
#' This may or may not be desirable, as it improves robustness to outliers but reduces the information available to obtain a highly resolved ranking. 
#' The most appropriate choice of effect size is left at the user's discretion.
#' 
#' Finally, the \code{logFC.detected} columns contain the log-fold change in the proportion of cells with detected (i.e., non-zero) expression between groups.
#' This is specifically useful for detecting binary expression patterns, e.g., activation of an otherwise silent gene.
#' Note that the non-zero status of the data is not affected by normalization, so differences in library size will implicitly affect the value of this metric.
#' However, this is not necessarily problematic for marker gene detection - users can treat this as \emph{retaining} information about the total RNA content, analogous to spike-in normalization.
#'
#' @section Computing effect size summaries:
#' To simplify interpretation, we summarize the effect sizes across all pairwise comparisons into a few key metrics.
#' For each group \eqn{X}, we consider the effect sizes from all pairwise comparisons between \eqn{X} and other groups. 
#' We then compute the following values:
#' \itemize{
#' \item \code{mean.*}, the mean effect sze across all pairwise comparisons involving \eqn{X}.
#' A large value (>0 for log-fold changes, >0.5 for the AUCs) indicates that the gene is upregulated in \eqn{X} compared to the average of the other groups.
#' \item \code{median.*}, the median effect size across all pairwise comparisons involving \eqn{X}.
#' A large value indicates that the gene is upregulated in \eqn{X} compared to most (>50\%) other groups.
#' \item \code{min.*}, the minimum effect size across all pairwise comparisons involving \eqn{X}.
#' A large value indicates that the gene is upregulated in \eqn{X} compared to all other groups,
#' while a small value (<0 for log-fold changes, <0.5 for the AUCs) indicates that the gene is downregulated in \eqn{X} compared to at least one other group.
#' \item \code{max.*}, the minimum effect size across all pairwise comparisons involving \eqn{X}.
#' A large value indicates that the gene is upregulated in \eqn{X} compared to at least one other group,
#' while a small value indicates that the gene is downregulated in \eqn{X} compared to all other groups.
#' }
#'
#' One set of these columns is added to the DataFrame for each effect size described above.
#' For example, the mean column for the AUC would be \code{mean.AUC}.
#' We can then reorder each group's DataFrame by our column of choice, depending on which summary and effect size we are interested in.
#' For example, if we ranked by decreasing \code{min.logFC.detected}, we would be aiming for marker genes that exhibit strong binary increases in expression in \eqn{X} compared to \emph{all} other groups.
#' 
#' The mean is obtained by weighting each comparison according to the number of cells involved.
#' This ensures that statistics from comparisons with very few cells do not skew the summary.
#' If \code{weight.fun=NULL}, the default weight is defined as the number of cells, capped at 100 cells;
#' this upweights comparisons with "enough" cells while ensuring that comparisons involving very large groups do not dominate the summary.
#' 
#' If \code{full.stats=TRUE}, an extra \code{full.*} column is returned in the DataFrame.
#' This contains a nested DataFrame with number of columns equal to the number of other groups.
#' Each column contains the statistic from the comparison between \eqn{X} and the other group.
#'
#' @author Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay, only using k-means for convenience.
#' kout <- kmeans(t(logcounts(sce)), centers=4) 
#'
#' out <- scoreMarkers(sce, groups=kout$cluster)
#' out
#'
#' of.interest <- out[[1]]
#' of.interest[order(of.interest$min.logFC.detected, decreasing=TRUE),]
#' 
#' @name scoreMarkers
NULL

#' @importFrom S4Vectors I DataFrame List
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom Matrix t
.scoreMarkers <- function(x, groups, block=NULL, design=NULL, weight.fun=NULL, row.data=NULL, full.stats=FALSE, BPPARAM=SerialParam()) {
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (!is.null(design)) {
        if (!is.null(block)) {
            stop("'block' and 'design' cannot both be specified")
        }

        # An unprincipled hack to deal with the design matrix for AUCs and the
        # logFC.detected. Especially for the latter; negative corrected values
        # are just set to zero and considered to be "undetected". Close enough.
        gdesign <- model.matrix(~0 + factor(groups))
        x <- ResidualMatrix::ResidualMatrix(
            t(DelayedArray(x)), 
            design=cbind(gdesign, design), 
            keep=seq_len(ncol(gdesign))
        )
        x <- t(x)
    }

    ugroups <- unique(groups)
    desired <- DataFrame(expand.grid(left=ugroups, right=ugroups))
    desired <- desired[desired$left != desired$right,,drop=FALSE]

    # Preparing the data structures for processing.
    combo.out <- .find_block_combinations(groups, block)
    combo.id <- combo.out$assignments
    combos <- combo.out$combinations

    comp.out <- .create_blocked_comparisons(combos, desired)
    left <- comp.out$left
    right <- comp.out$right
    comparisons <- comp.out$comparisons

    ncells <- tabulate(combo.id, nbins=nrow(combos))
    left.n <- ncells[left]
    right.n <- ncells[right]

    involved <- seq_len(nrow(combos)) %in% left | seq_len(nrow(combos)) %in% right

    # Performing the per-cell calculations and gathering the statistics.
    stats <- rowBlockApply(x, combo.id=combo.id, 
        left=left, right=right, left.n=left.n, right.n=right.n, involved=involved,
        BPPARAM=BPPARAM)

    res <- do.call(mapply, c(list(FUN=rbind), stats, list(SIMPLIFY=FALSE, USE.NAMES=FALSE)))
    names(res) <- names(stats[[1]])

    # Averaging across blocks and then collating.
    merge.combos <- .create_block_mergers(combos)
    w <- left.n * right.n
    for (i in seq_along(res)) {
        effect.name <- names(res)[i]
        if (effect.name=="AUC") {
            REVERSE <- function(x) 1-x
        } else {
            REVERSE <- function(x) -x
        }

        ave.out <- .combine_across_blocks(res[[i]], weights=w, to.merge=merge.combos)

        res[[i]] <- .collate_into_DataFrame(comparisons, 
            averaged=ave.out$averaged, 
            weights=ave.out$weights, 
            desired=desired, 
            REVERSE=REVERSE,
            effect.name=effect.name,
            row.names=rownames(x), 
            full=full.stats) 
    }
    
    do.call(cbind, res)
}

#####################################################################
#####################################################################

#' @importFrom S4Vectors selfmatch
.uniquify_DataFrame <- function(df) {
    id <- selfmatch(df)
    keep <- !duplicated(id)
    uid <- id[keep]
    ucombos <- df[keep,,drop=FALSE]
    list(unique=ucombos, id=match(id, uid))
}

.find_block_combinations <- function(group, block) {
    everything <- DataFrame(group=groups)
    if (!is.null(block)) {
        everything$block <- block
    }
    out <- .uniquify_DataFrame(everything)
    list(combinations=out$unique, assignments=out$id)
}

#' @importFrom S4Vectors DataFrame
.create_blocked_comparisons <- function(ucombos, desired) {
    rows <- seq_len(nrow(ucombos))
    if (!is.null(ucombos$block)) {
        by.block <- split(rows, ucombos$block)
    } else {
        by.block <- list(rows)
    }

    left <- right <- actual.comp <- vector("list", length(by.block))
    for (i in seq_along(by.block)) {
        current <- by.block[[i]]

        all.pairs <- DataFrame(expand.grid(left=current, right=current))
        all.pairs <- all.pairs[all.pairs %in% desired,,drop=FALSE]
        keep <- all.pairs$left > all.pairs$right
        all.pairs <- all.pairs[keep,]

        left[[i]] <- all.pairs$left
        right[[i]] <- all.pairs$right
    }

    block.id <- rep(seq_along(left), lengths(left))
    if (!is.null(names(by.block))) {
        block.id <- names(by.block)[block.id]
    }

    left <- unlist(left)
    right <- unlist(right)
    actual.comp <- DataFrame(left=ucombos$group[left], right=ucombos$group[right], block=block.id)

    list(left=left, right=right, comparisons=actual.comp)
}

#####################################################################
#####################################################################

#' @importFrom scuttle summarizeAssayByGroup
.compute_all_effect_sizes <- function(x, combo.id, left, right, left.n, right.n, involved) {
    stats <- .compute_mean_var(x, BPPARAM=SerialParam(), subset.row=NULL, design=NULL, block.FUN=compute_blocked_stats_none, block=combo.id)
    cohen <- .compute_pairwise_cohen_d(stats$means, stats$vars, stats$ncells, left, right, left.n, right.n)

    detected <- summarizeAssayByGroup(x, combo.id, statistics="num.detected")
    o <- order(detected$id)
    lfc <- .compute_lfc_detected(assay(detected, withDimnames=FALSE)[,o,drop=FALSE],
        detected$ncells[o], left, right, left.n, right.n)

    auc <- .compute_auc(x, combo.id, involved, left, right, left.n, right.n)

    list(logFC.cohen=cohen, AUC=auc, logFC.detected=lfc)
}

#' @importFrom DelayedMatrixStats rowWeightedMeans
.compute_pairwise_cohen_d <- function(means, vars, ncells, left, right, left.n, right.n) {
    all.delta <- means[,left,drop=FALSE] - stats$means[,right,drop=FALSE]

    left.s2 <- t(vars[,left,drop=FALSE])
    right.s2 <- t(vars[,right,drop=FALSE])
    pooled.s2 <- ((left.n - 1) * left.s2 + (right.n - 1) * right.s2)/(left.n + right.n - 2)
    pooled.s2 <- t(pooled.s2)

    is.zero <- all.delta == 0
    d <- all.delta / sqrt(pooled.s2)
    d[is.zero] <- 0

    d
}

.compute_lfc_detected <- function(ndetected, ncells, left, right, left.n, right.n) {
    left.detected <- ndetected[,left,drop=FALSE]
    right.detected <- ndetected[,right,drop=FALSE]

    mean.n <- (left.n + right.n)/2
    left.pseudo <- 1 * left.n / mean.n
    right.pseudo <- 1 * right.n / mean.n

    left.prop <- (t(left.detected) + left.pseudo) / (left.n + left.pseudo * 2)
    right.prop <- (t(right.detected) + right.pseudo) / (right.n + right.pseudo * 2)
    log2(t(left.prop/right.prop))
}

.compute_auc <- function(x, combinations, involved, left, right, left.n, right.n) {
    overlap <- overlap_exprs_paired(x, combinations, left, right, involved)
    t(t(overlap) / (left.n * right.n))
}

#####################################################################
#####################################################################

#' @importFrom S4Vectors splitAsList
.create_block_mergers <- function(combinations) {
    u.out <- .uniquify_DataFrame(combinations[,1:2])
    f <- factor(u.out$id, seq_len(nrow(u.out$unique)))
    list(comparisons=u.out$unique, block=splitAsList(seq_along(f), f))
}

.combine_across_blocks <- function(statistics, weights, to.merge) {
    output <- vector("list", length(to.merge))
    combined.weight <- numeric(length(to.merge))

    for (i in seq_along(to.merge)) {
        current <- statistics[,to.merge[[i]],drop=FALSE]
        ns <- sum(weights)
        combined.weight[[i]] <- ns
        output[[i]] <- colSums(t(current) * weights, na.rm=TRUE) / ns
    }

    list(averaged=output, weights=combined.weight)
}

#' @importFrom S4Vectors splitAsList DataFrame
#' @importFrom DelayedMatrixStats rowMins rowMedians rowWeightedMeans
.collate_into_DataFrame <- function(comparisons, averaged, weights, desired, REVERSE, effect.name, row.names=NULL, full=FALSE) {
    m <- match(desired, comparisons)

    flipped <- comparisons[,2:1]
    colnames(flipped) <- colnames(comparisons)
    fm <- match(desired, flipped)

    by.left <- splitAsList(seq_len(nrow(desired)), desired$left)
    for (i in names(by.left)) {
        chosen <- by.left[[i]]
        collated <- vector("list", length(chosen))
        names(collated) <- desired$right[chosen]
        w <- numeric(length(chosen))

        original <- m[chosen]
        keep <- !is.na(original)
        original.kept <- original[keep]
        collated[keep] <- averaged[original.kept]
        w[keep] <- weights[original.kept]

        flipped <- fm[chosen]
        fkeep <- !is.na(flipped)
        flipped.kept <- flipped[keep]
        collated[fkeep] <- lapply(averaged[flipped.kept], REVERSE)
        w[fkeep] <- weights[flipped.kept]

        mat <- do.call(cbind, collated)
        df <- DataFrame(
            mean=rowWeightedMeans(effect.mat, w=w, na.rm=TRUE),
            min=rowMins(effect.mat, na.rm=TRUE),
            median=rowMedians(effect.mat, na.rm=TRUE),
            max=rowMaxs(effect.mat, na.rm=TRUE),
            row.names=rn
        )
        if (full) {
            df$full <- DataFrame(collated, row.names=rn)
        }
        colnames(df) <- paste0(colnames(df), ".", effect.name)

        by.left[[i]] <- df
    }

    by.left
}

#####################################################################
#####################################################################

#' @export
#' @rdname scoreMarkers
setGeneric("scoreMarkers", function(x, groups, ...) standardGeneric("scoreMarkers"))

#' @export
#' @rdname scoreMarkers
setMethod("scoreMarkers", "ANY", .scoreMarkers) 

#' @export
#' @rdname scoreMarkers
#' @importFrom SummarizedExperiment assay
setMethod("scoreMarkers", "SummarizedExperiment", function(x, groups, ..., assay.type="logcounts") {
    .scoreMarkers(assay(x, assay.type), groups, ...)
})

#' @export
#' @rdname scoreMarkers
#' @importFrom SingleCellExperiment colLabels
setMethod("scoreMarkers", "SingleCellExperiment", function(x, groups=colLabels(x, onAbsence="error"), ...) {
    callNextMethod()
})
