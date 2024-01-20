#' Score marker genes
#'
#' Compute various summary scores for potential marker genes to distinguish between groups of cells.
#'
#' @param x A matrix-like object containing log-normalized expression values, with genes in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix in its assays.
#' @param groups A factor or vector containing the identity of the group for each cell in \code{x}.
#' @param block A factor or vector specifying the blocking level for each cell in \code{x}.
#' @param row.data A DataFrame with the same number and names of rows in \code{x}, containing extra information to insert into each DataFrame.
#' @param full.stats Logical scalar indicating whether the statistics from the pairwise comparisons should be directly returned.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the calculations should be parallelized.
#' @param ... For the generic, further arguments to pass to individual methods.
#' 
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param assay.type String or integer scalar specifying the assay containing the log-expression matrix to use.
#' @param lfc Numeric scalar specifying the log-fold change threshold to compute effect sizes against.
#' @param pairings A vector, list or matrix specifying how the comparisons are to be performed, see details.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
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
#' @section Choice of effect sizes:
#' The \code{logFC.cohen} columns contain the standardized log-fold change, i.e., Cohen's d.
#' For each pairwise comparison, this is defined as the difference in the mean log-expression for each group scaled by the average standard deviation across the two groups.
#' (Technically, we should use the pooled variance; however, this introduces some unpleasant asymmetry depending on the variance of the larger group, so we take a simple average instead.)
#' Cohen's d is analogous to the t-statistic in a two-sample t-test and avoids spuriously large effect sizes from comparisons between highly variable groups.
#' We can also interpret Cohen's d as the number of standard deviations between the two group means.
#' 
#' The \code{AUC} columns contain the area under the curve.
#' This is the probability that a randomly chosen observation in one group is greater than a randomly chosen observation in the other group.
#' The AUC is closely related to the U-statistic used in the Wilcoxon rank sum test.
#' Values greater than 0.5 indicate that a gene is upregulated in the first group.
#'
#' The key difference between the AUC and Cohen's d is that the former is less sensitive to the variance within each group.
#' The clearest example is that of two distributions that exhibit no overlap, where the AUC is the same regardless of the variance of each distribution.
#' This may or may not be desirable as it improves robustness to outliers but reduces the information available to obtain a highly resolved ranking. 
#' The most appropriate choice of effect size is left at the user's discretion.
#'
#' Finally, the \code{logFC.detected} columns contain the log-fold change in the proportion of cells with detected (i.e., non-zero) expression between groups.
#' This is specifically useful for detecting binary expression patterns, e.g., activation of an otherwise silent gene.
#' Note that the non-zero status of the data is not altered by normalization, so differences in library size will not be removed when computing this metric.
#' This effect is not necessarily problematic - users can interpret it as \emph{retaining} information about the total RNA content, analogous to spike-in normalization.
#'
#' @section Setting a log-fold change threshold:
#' The default settings may yield highly ranked genes with large effect sizes but low log-fold changes if the variance is low (Cohen's d) or separation is clear (AUC).
#' Such genes may not be particularly interesting as the actual change in expression is modest.
#' Setting \code{lfc} allows us to focus on genes with large log-fold changes between groups,
#' by simply shifting the \dQuote{other} group's expression values by \code{lfc} before computing effect sizes.
#'
#' When \code{lfc} is not zero, Cohen's d is generalized to the standardized difference between the observed log-fold change and \code{lfc}.
#' For example, if we had \code{lfc=2} and we obtained a Cohen's d of 3, this means that the observed log-fold change was 3 standard deviations above a value of 2.
#' A side effect is that we can only unambiguously interpret the direction of Cohen's d when it has the same sign as \code{lfc}.
#' Our above example represents upregulation, but if our Cohen's d was negative, this could either mean downregulation or simply that our observed log-fold change was less than \code{lfc}.
#' 
#' When \code{lfc} is not zero, the AUC is generalized to the probability of obtaining a random observation in one group that is greater than a random observation plus \code{lfc} in the other group.
#' For example, if we had \code{lfc=2} and we obtained an AUC of 0.8, this means that we would observe a difference of \code{lfc} or greater between the random observations.
#' Again, we can only unambiguously interpret the direction of the change when it is the same as the sign of the \code{lfc}.
#' In this case, an AUC above 0.5 with a positive \code{lfc} represents upregulation, but an AUC below 0.5 could mean either downregulation or a log-fold change less than \code{lfc}.
#' 
#' A non-zero setting of \code{lfc} has no effect on the log-fold change in the proportion of cells with detected expression.
#' 
#' @section Computing effect size summaries:
#' To simplify interpretation, we summarize the effect sizes across all pairwise comparisons into a few key metrics.
#' For each group \eqn{X}, we consider the effect sizes from all pairwise comparisons between \eqn{X} and other groups. 
#' We then compute the following values:
#' \itemize{
#' \item \code{mean.*}, the mean effect size across all pairwise comparisons involving \eqn{X}.
#' A large value (>0 for log-fold changes, >0.5 for the AUCs) indicates that the gene is upregulated in \eqn{X} compared to the average of the other groups.
#' A small value (<0 for the log-fold changes, <0.5 for the AUCs) indicates that the gene is downregulated in \eqn{X} instead.
#' \item \code{median.*}, the median effect size across all pairwise comparisons involving \eqn{X}.
#' A large value indicates that the gene is upregulated in \eqn{X} compared to most (>50\%) other groups.
#' A small value indicates that the gene is downregulated in \eqn{X} instead.
#' \item \code{min.*}, the minimum effect size across all pairwise comparisons involving \eqn{X}.
#' A large value indicates that the gene is upregulated in \eqn{X} compared to all other groups.
#' A small value indicates that the gene is downregulated in \eqn{X} compared to at least one other group.
#' \item \code{max.*}, the maximum effect size across all pairwise comparisons involving \eqn{X}.
#' A large value indicates that the gene is upregulated in \eqn{X} compared to at least one other group.
#' A small value indicates that the gene is downregulated in \eqn{X} compared to all other groups.
#' \item \code{rank.*}, the minimum rank (i.e., \dQuote{min-rank}) across all pairwise comparisons involving \eqn{X} - see \code{?\link{computeMinRank}} for details.
#' A small min-rank indicates that the gene is one of the top upregulated genes in at least one comparison to another group.
#' }
#' One set of these columns is added to the DataFrame for each effect size described above.
#' For example, the mean column for the AUC would be \code{mean.AUC}.
#' We can then reorder each group's DataFrame by our column of choice, depending on which summary and effect size we are interested in.
#' For example, if we ranked by decreasing \code{min.logFC.detected}, we would be aiming for marker genes that exhibit strong binary increases in expression in \eqn{X} compared to \emph{all} other groups.
#' 
#' If \code{full.stats=TRUE}, an extra \code{full.*} column is returned in the DataFrame.
#' This contains a nested DataFrame with number of columns equal to the number of other groups.
#' Each column contains the statistic from the comparison between \eqn{X} and the other group.
#'
#' Keep in mind that the interpretations above also depend on the sign of \code{lfc}.
#' The concept of a \dQuote{large} summary statistic (>0 for Cohen's d, >0.5 for the AUCs) can only be interpreted as upregulation when \code{lfc >= 0}.
#' Similarly, the concept of a \dQuote{small} value (<0 for Cohen's d, <0.5 for the AUCs) cannot be interpreted as downregulation when \code{lfc <= 0}.
#' For example, if \code{lfc=1}, a positive \code{min.logFC.cohen} can still be interpreted as upregulation in \eqn{X} compared to all other groups,
#' but a negative \code{max.logFC.cohen} could not be interpreted as downregulation in \eqn{X} compared to all other groups.
#'
#' @section Computing other descriptive statistics: 
#' We report the mean log-expression of all cells in \eqn{X}, as well as the grand mean of mean log-expression values for all other groups.
#' This is purely descriptive; while it can be used to compute an overall log-fold change, ranking is best performed on one of the effect sizes described above.
#' We also report the proportion of cells with detectable expression in \eqn{X} and the mean proportion for all other groups.
#'
#' When \code{block} is specified, the reported mean for each group is computed via \code{\link{correctGroupSummary}}. 
#' Briefly, this involves fitting a linear model to remove the effect of the blocking factor from the per-group mean log-expression.
#' The same is done for the detected proportion, except that the values are subjected to a logit transformation prior to the model fitting.
#' In both cases, each group/block combination is weighted by its number of cells in the model.
#'
#' @section Controlling the pairings:
#' The \code{pairings} argument specifies the pairs of groups that should be compared.
#' This can be:
#' \itemize{
#' \item \code{NULL}, in which case comparisons are performed between all groups in \code{groups}.
#' \item A vector of the same type as \code{group}, specifying a subset of groups of interest.
#' We then perform all pairwise comparisons between groups in the subset.
#' \item A list of two vectors, each of the same type as \code{group} and specifying a subset of groups.
#' Comparisons are performed between one group from the first vector and another group from the second vector.
#' \item A matrix of two columns of the same type as \code{group}.
#' Each row is assumed to specify a pair of groups to be compared.
#' }
#' 
#' Effect sizes (and their summaries) are computed for only the pairwise comparisons specified by \code{pairings}.
#' Similarly, the \code{other.*} values in \eqn{X}'s DataFrame are computed using only the groups involved in pairwise comparisons with \eqn{X}.
#' The default of \code{pairings=NULL} ensures that all groups are used and effect sizes for all pairwise comparisons are reported;
#' however, this may not be the case for other values of \code{pairings}.
#' 
#' For list and matrix arguments, the first vector/column is treated as the first group in the effect size calculations.
#' Statistics for each comparison will only show up in the DataFrame for the first group, 
#' i.e., a comparison between \eqn{X} and \eqn{Y} will have a valid \code{full.AUC$Y} field in \eqn{X}'s DataFrame but not vice versa.
#' If both directions are desired in the output, both of the corresponding permutations should be explicitly specified in \code{pairings}.
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
#' # Ranking by a metric of choice:
#' of.interest <- out[[1]]
#' of.interest[order(of.interest$mean.AUC, decreasing=TRUE),1:4]
#' of.interest[order(of.interest$rank.AUC),1:4]
#' of.interest[order(of.interest$median.logFC.cohen, decreasing=TRUE),1:4]
#' of.interest[order(of.interest$min.logFC.detected, decreasing=TRUE),1:4]
#' 
#' @name scoreMarkers
NULL

#' @importFrom S4Vectors I DataFrame SimpleList
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel SerialParam 
#' @importFrom Matrix t
.scoreMarkers <- function(x, groups, block=NULL, pairings=NULL, lfc=0, row.data=NULL, full.stats=FALSE, subset.row=NULL, BPPARAM=SerialParam()) {
#    if (!is.null(design)) {
#        if (!is.null(block)) {
#            stop("'block' and 'design' cannot both be specified")
#        }
#
#        # An unprincipled hack to deal with the design matrix for AUCs and the
#        # logFC.detected. Especially for the latter; negative corrected values
#        # are just set to zero and considered to be "undetected". Close enough.
#        gdesign <- model.matrix(~0 + factor(groups))
#        x <- ResidualMatrix::ResidualMatrix(
#            t(DelayedArray(x)), 
#            design=cbind(gdesign, design), 
#            keep=seq_len(ncol(gdesign))
#        )
#        x <- t(x)
#    }

    # Define the desired comparisons here.
    universe <- sort(unique(groups))
    desired.out <- .expand_pairings(pairings, universe=universe)
    desired.comparisons <- DataFrame(left=universe[desired.out$id1], right=universe[desired.out$id2])

    keep <- groups %in% desired.comparisons[,1] | groups %in% desired.comparisons[,2]
    if (!all(keep)) {
        # Trimming the matrix to avoid redundant calculations.
        x <- x[,keep,drop=FALSE]
        groups <- groups[keep]
        block <- block[keep]
    }

    # Preparing the data structures for processing.
    combination.out <- .group_block_combinations(groups, block)
    combination.id <- combination.out$id
    unique.combinations <- combination.out$combinations
    ncells <- as.double(tabulate(combination.id, nbins=nrow(unique.combinations)))

    collapse.symmetric <- lfc==0
    reindexed.comparisons <- .reindex_comparisons_for_combinations(unique.combinations, desired.comparisons, collapse.symmetric=collapse.symmetric)
    left <- reindexed.comparisons$left
    right <- reindexed.comparisons$right

    involved <- .group_by_used_combinations(combination.id, left, right, nrow(unique.combinations))
    pre.ave <- .identify_effects_to_average(unique.combinations, reindexed.comparisons)
    desired.indices <- .cross_reference_to_desired(pre.ave$averaged.comparisons, desired.comparisons, collapse.symmetric=collapse.symmetric)

    # Cleaning up the rows.
    if (!is.null(row.data)) {
        if (!identical(nrow(row.data), nrow(x))) {
            stop("'row.data' and 'x' should have the same number of rows")
        }
        if (!identical(rownames(row.data), rownames(x))) {
            stop("'row.data' and 'x' should have the same row names")
        }
    }
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
        row.data <- row.data[subset.row,,drop=FALSE]
    }

    # Performing the per-cell calculations and gathering the statistics.
    stats <- rowBlockApply(x, 
        FUN=.compute_all_effect_sizes,
        combination.id=combination.id, 
        left=left, 
        right=right, 
        ncells=ncells,
        unique.combinations=unique.combinations, 
        indices.to.average=pre.ave$indices.to.average, 
        desired.indices=desired.indices,
        involved=involved,
        lfc=lfc,
        full.stats=full.stats,
        BPPARAM=BPPARAM)

    res <- .mapply_bind(stats, rbind)

    # Computing the minimum rank across genes.
    for (i in names(res)) {
        all.names <- grep("^rank\\.", colnames(res[[i]]))
        for (j in all.names) {
            res[[i]][[j]] <- computeMinRank(res[[i]][[j]])
        }
    }

    # Slapping on the row data.
    if (!is.null(row.data)) {
        for (i in seq_along(res)) {
            res[[i]] <- cbind(row.data, res[[i]])
        }
    }

    SimpleList(res)
}

#####################################################################
#####################################################################

#' @importFrom S4Vectors selfmatch
.uniquify_DataFrame <- function(df) {
    id <- selfmatch(df)
    keep <- !duplicated(id)
    uid <- id[keep]

    # Ordering to ensure that earlier combination IDs imply earlier DF values.
    # This ensures that we can deduplicate comparisons by removing those where
    # the left ID < right ID, which implies that the left group < right group.
    # Otherwise, if the combinations were not ordered, we would have some 
    # comparisons where the left group < right group and others where the 
    # left group > right group, despite all of them having left ID < right ID;
    # this prevents collapsing of the same comparison across multiple blocks.
    ucombos <- df[keep,,drop=FALSE]
    o <- order(ucombos)
    uid <- uid[o]
    ucombos <- ucombos[o,,drop=FALSE]

    list(unique=ucombos, id=match(id, uid))
}

#' @importFrom S4Vectors DataFrame
.group_block_combinations <- function(groups, block) {
    everything <- DataFrame(group=groups)
    if (!is.null(block)) {
        everything$block <- block
    }
    out <- .uniquify_DataFrame(everything)
    list(combinations=out$unique, id=out$id)
}

#' @importFrom S4Vectors DataFrame
#' @importMethodsFrom S4Vectors %in%
.reindex_comparisons_for_combinations <- function(unique.combinations, desired.comparisons, collapse.symmetric=TRUE) {
    rows <- seq_len(nrow(unique.combinations))
    if (!is.null(unique.combinations$block)) {
        by.block <- split(rows, unique.combinations$block)
    } else {
        by.block <- list(rows)
    }

    left <- right <- block.val <- vector("list", length(by.block))
    for (i in seq_along(by.block)) {
        indices <- by.block[[i]]
        group <- unique.combinations$group[indices]

        # Rephrasing the desired comparisons in terms of the groups available in this block.
        mleft <- match(desired.comparisons$left, group)
        mright <- match(desired.comparisons$right, group)
        keep <- !is.na(mleft) & !is.na(mright)

        index.left <- indices[mleft[keep]]
        index.right <- indices[mright[keep]]

        if (collapse.symmetric) { 
            # Eliminating redundant comparisons with flipped orientations. This
            # is only permissible when the effect sizes are somehow symmetric.
            index.pairs <- DataFrame(left=pmax(index.left, index.right), right=pmin(index.left, index.right))
        } else {
            index.pairs <- DataFrame(left=index.left, right=index.right)
        }
        index.pairs <- unique(index.pairs)

        left[[i]] <- index.pairs$left
        right[[i]] <- index.pairs$right
    }

    DataFrame(left=unlist(left), right=unlist(right))
}

.group_by_used_combinations <- function(combination.id, left, right, ncombinations) {
    useful <- which(combination.id %in% c(left, right))
    f <- factor(combination.id[useful], seq_len(ncombinations))
    split(useful - 1L, f)
}

#####################################################################
#####################################################################

#' @importFrom scuttle summarizeAssayByGroup correctGroupSummary
.compute_all_effect_sizes <- function(x, combination.id, left, right, ncells,
    involved, unique.combinations, indices.to.average, desired.indices, lfc, full.stats)
{
    left.ncells <- ncells[left]
    right.ncells <- ncells[right]

    # Computing effects.
    stats <- compute_blocked_stats_none(x, combination.id - 1L, length(involved))
    cohen <- .compute_pairwise_cohen_d(stats[[1]], stats[[2]], left, right, lfc=lfc)

    detected.se <- summarizeAssayByGroup(x, combination.id, statistics=c("num.detected", "prop.detected"))
    m <- match(seq_len(nrow(unique.combinations)), detected.se$ids)
    ndetected <- assay(detected.se, withDimnames=FALSE)[,m,drop=FALSE]
    nlfc <- .compute_lfc_detected(ndetected, left, right, left.ncells, right.ncells)

    auc <- .compute_auc(x, involved, left, right, left.ncells, right.ncells, lfc=lfc)

    # Averaging across blocks and then collating.
    weights <- left.ncells * right.ncells
    output <- list(logFC.cohen=cohen, AUC=auc, logFC.detected=nlfc)

    for (effect in names(output)) {
        if (effect == "AUC") {
            REVERSE <- function(x) 1-x
        } else {
            REVERSE <- function(x) -x
        }

        ave.out <- .average_effect_across_blocks(output[[effect]], weights=weights, indices.to.average=indices.to.average)

        output[[effect]] <- .collate_into_DataFrame(
            desired.indices,
            averaged.effects=ave.out$averaged.effects, 
            combined.weights=ave.out$combined.weights, 
            REVERSE=REVERSE,
            effect.name=effect,
            nrow=nrow(x),
            row.names=rownames(x), 
            full.stats=full.stats) 
    }

    output <- .mapply_bind(output, cbind)

    # Computing mean stats across blocks.
    FUN <- function(...) correctGroupSummary(..., group=unique.combinations$group, block=unique.combinations$block, weights=ncells)
    more.stats <- list(
        detected=FUN(assay(detected.se, "prop.detected", withDimnames=FALSE)[,m,drop=FALSE], transform="logit"),
        average=FUN(x=stats[[1]])
    )

    for (i in names(more.stats)) {
        mat <- more.stats[[i]]
        for (j in names(output)) {
            existing <- output[[j]]
            self <- colnames(mat) == j
            other <- colnames(mat) %in% as.character(desired.indices[[j]]$right)
            extra <- DataFrame(self=rowMeans(mat[,self,drop=FALSE]), other=rowMeans(mat[,other,drop=FALSE]), row.names=rownames(existing))
            colnames(extra) <- paste0(colnames(extra), ".", i)
            output[[j]] <- cbind(extra, existing)
        }
    }

    output
}

#' @importFrom DelayedMatrixStats rowWeightedMeans
.compute_pairwise_cohen_d <- function(means, vars, left, right, lfc) {
    all.delta <- means[,left,drop=FALSE] - means[,right,drop=FALSE] - lfc
    pooled.s2 <- (vars[,left,drop=FALSE] + vars[,right,drop=FALSE])/2

    is.zero <- all.delta == 0
    d <- all.delta / sqrt(pooled.s2)
    d[is.zero] <- 0

    d
}

.compute_lfc_detected <- function(ndetected, left, right, left.ncells, right.ncells) {
    left.detected <- ndetected[,left,drop=FALSE]
    right.detected <- ndetected[,right,drop=FALSE]

    # Using the minimum to provide greater shrinkage when either group has very few cells.
    min.n <- pmin(left.ncells, right.ncells)
    left.pseudo <- 1 * left.ncells / min.n
    right.pseudo <- 1 * right.ncells / min.n

    left.prop <- (t(left.detected) + left.pseudo) / (left.ncells + left.pseudo * 2)
    right.prop <- (t(right.detected) + right.pseudo) / (right.ncells + right.pseudo * 2)
    log2(t(left.prop/right.prop))
}

.compute_auc <- function(x, involved, left, right, left.ncells, right.ncells, lfc) {
    overlap <- overlap_exprs_paired(x, left, right, involved, lfc)
    t(overlap / (left.ncells * right.ncells))
}


#####################################################################
#####################################################################

.mapply_bind <- function(df.list, FUN) {
    ref <- df.list[[1]]
    for (i in seq_along(df.list)[-1]) {
        chosen <- df.list[[i]]
        for (j in seq_along(chosen)) {
            ref[[j]] <- FUN(ref[[j]], chosen[[j]])
        }
    }
    ref
}

#' @importFrom S4Vectors DataFrame splitAsList
.identify_effects_to_average <- function(unique.combinations, reindexed.comparisons) {
    df <- DataFrame(
        left = unique.combinations$group[reindexed.comparisons$left],
        right = unique.combinations$group[reindexed.comparisons$right]
    )

    u.out <- .uniquify_DataFrame(df)
    comp <- u.out$unique

    f <- factor(u.out$id, seq_len(nrow(comp)))
    merger <- splitAsList(seq_along(f), f)

    list(averaged.comparisons=comp, indices.to.average=merger)
}

#' @importFrom DelayedMatrixStats rowWeightedMeans 
.average_effect_across_blocks <- function(effects, weights, indices.to.average) {
    output <- vector("list", length(indices.to.average))
    combined.weight <- numeric(length(indices.to.average))

    for (i in seq_along(indices.to.average)) {
        m <- indices.to.average[[i]]
        cur.e <- effects[,m,drop=FALSE]
        cur.w <- weights[m]
        combined.weight[[i]] <- sum(cur.w)
        output[[i]] <- rowWeightedMeans(cur.e, cur.w, na.rm=TRUE)
    }

    list(averaged.effects=output, combined.weights=combined.weight)
}

#' @importFrom S4Vectors match
.cross_reference_to_desired <- function(averaged.comparisons, desired.comparisons, collapse.symmetric=TRUE) {
    m <- match(desired.comparisons, averaged.comparisons)
   
    if (collapse.symmetric) {
        # If the symmetric comparisons were deduplicated by flipping left/right,
        # we flip them for the cross-referencing back to 'desired.comparisons'.
        flipped.comparisons <- averaged.comparisons[,2:1]
        colnames(flipped.comparisons) <- colnames(averaged.comparisons)
        fm <- match(desired.comparisons, flipped.comparisons)
    } else {
        fm <- rep(NA_integer_, nrow(desired.comparisons))
    }

    # drop=TRUE to ignore unused factor levels, which are highly unlikely to be desired.
    indices <- split(seq_len(nrow(desired.comparisons)), desired.comparisons$left, drop=TRUE)

    for (i in seq_along(indices)) {
        chosen <- indices[[i]]
        indices[[i]] <- list(
            right=desired.comparisons$right[chosen],
            direct.match=m[chosen],
            flipped.match=fm[chosen]
        )
    }

    indices
}

#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom DelayedMatrixStats rowMins rowMedians rowMaxs 
#' @importFrom Matrix rowMeans
.collate_into_DataFrame <- function(desired.indices, averaged.effects, combined.weights, REVERSE, effect.name, nrow, row.names=NULL, full.stats=FALSE) {
    output <- desired.indices

    for (i in names(desired.indices)) {
        current <- desired.indices[[i]]

        right <- current$right
        collated <- vector("list", length(right))
        names(collated) <- right
        w <- numeric(length(right))

        original <- current$direct.match
        keep <- !is.na(original)
        if (any(keep)) {
            original.kept <- original[keep]
            collated[keep] <- averaged.effects[original.kept]
            w[keep] <- combined.weights[original.kept]
        }

        flipped <- current$flipped.match
        fkeep <- !is.na(flipped)
        if (any(fkeep)) {
            flipped.kept <- flipped[fkeep]
            collated[fkeep] <- lapply(averaged.effects[flipped.kept], REVERSE)
            w[fkeep] <- combined.weights[flipped.kept]
        }

        # Filling in the leftovers, e.g., due to an impossible comparison.
        for (j in which(!keep & !fkeep)) {
            collated[[j]] <- rep(NA_real_, nrow)
        }

        # 'current', and thus 'collated', is guaranteed to be non-empty due to drop=TRUE.
        # As such, there's no need to implement any protection in the constructor.
        full <- DataFrame(collated, row.names=row.names, check.names=FALSE)

        effect.mat <- as.matrix(full)
        df <- DataFrame(
            mean=rowMeans(effect.mat, na.rm=TRUE),
            min=rowMins(effect.mat, na.rm=TRUE),
            median=rowMedians(effect.mat, na.rm=TRUE),
            max=rowMaxs(effect.mat, na.rm=TRUE),
            rank=I(effect.mat), # need to combine the full effect size matrix across rows before compute ranks.
            row.names=row.names
        )

        if (full.stats) {
            mcols(full)$weight <- w
            df$full <- full
        }

        colnames(df) <- paste0(colnames(df), ".", effect.name)

        output[[i]] <- df
    }

    output
}

#####################################################################
#####################################################################

#' @export
#' @rdname scoreMarkers
setGeneric("scoreMarkers", function(x, ...) standardGeneric("scoreMarkers"))

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
    callNextMethod(x, groups=groups, ...)
})
