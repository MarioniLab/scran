#' Score marker genes
#'
#' Compute various summary scores for potential marker genes to distinguish between groups of cells.
#'
#' @param x A matrix-like object containing log-normalized expression values, with genes in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix in its assays.
#' @param groups A factor or vector containing the identity of the group for each cell in \code{x}.
#' @param ... Further arguments to be passed to \code{\link{pairwiseTTest}} and related functions.
#' The most obvious of these is \code{block}.
#' @param weight.fun Function indicating how the statistics from different comparisons should be weighted.
#' This should accept a vector of integers containing the number of cells involved in the comparison,
#' and return a vector of equal length containing the weights. 
#' Defaults to \code{log10}.
#' @param full.stats Logical scalar indicating whether the statistics from the pairwise comparisons should be directly returned.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the calculations should be parallelized.
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
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to the average of the other groups.
#' \item \code{median.*}, the median effect size across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to most (>50\%) other groups.
#' \item \code{min.*}, the minimum effect size across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to all other groups,
#' while a negative value indicates that the gene is downregulated in \eqn{X} compared to at least one other group.
#' \item \code{max.*}, the minimum effect size across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to at least one other group,
#' while a negative value indicates that the gene is downregulated in \eqn{X} compared to all other groups.
#' }
#'
#' One set of these columns is added to the DataFrame for each effect size described above.
#' For example, the mean column for the AUC would be \code{mean.AUC}.
#' We can then reorder each group's DataFrame by our column of choice, depending on which summary and effect size we are interested in.
#' For example, if we ranked by decreasing \code{min.logFC.detected}, we would be aiming for marker genes that exhibit strong binary increases in expression in \eqn{X} compared to \emph{all} other groups.
#' 
#' In practice, the mean and the median metrics are obtained by weighting each comparison according to the number of cells involved.
#' This ensures that statistics from comparisons with very few cells do not skew the summaries.
#' If \code{weight.fun=NULL}, the default weight is defined as the log10-transformed number of cells, capped at 100 cells;
#' this favors comparisons with "enough" cells while ensuring that comparisons involving very large groups do not dominate the summary.
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
#' @export
#' @importFrom S4Vectors I DataFrame List
#' @importFrom MatrixGenerics rowMins rowWeightedMedians rowWeightedMeans
scoreMarkers <- function(x, groups, ..., weight.fun=NULL, full.stats=FALSE, BPPARAM=SerialParam()) {
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # TODO: avoid the need to do multiple extractions from 'x' 
    # when computing gene-wise statistics.
    output.p <- output.effect <- list()

    for (tt in c("t", "wilcox", "binom")) {
        if (tt == "t") {
            FUN <- function(...) pairwiseTTests(..., std.lfc=TRUE)
            effect.field <- "logFC"
            output.field <- "logFC.cohen"
        } else if (tt == "wilcox") {
            FUN <- pairwiseWilcox
            effect.field <- output.field <- "AUC"
        } else {
            FUN <- pairwiseBinom
            effect.field <- "logFC"
            output.field <- "logFC.detected"
        }

        fit <- FUN(x, groups, ..., log.p=TRUE, BPPARAM=BPPARAM)
        output <- combineMarkers(fit$statistics, fit$pairs, 
            log.p.in=TRUE, log.p.out=TRUE, full.stats=TRUE, pval.field="log.p.value", 
            effect.field=effect.field, sorted=FALSE, BPPARAM=BPPARAM)

        current.p <- current.effect <- list()
        for (i in names(output)) {
            df <- output[[i]]
            my.stats <- as.list(df[,grepl("^stats.", colnames(df)),drop=FALSE])
            my.effects <- lapply(my.stats, function(x) x[[effect.field]])
            names(my.effects) <- sub("^stats.", "", names(my.effects))
            current.effect[[i]] <- DataFrame(my.effects, row.names=rownames(df), check.names=FALSE)
        }

        output.effect[[output.field]] <- current.effect
    }

    # Computing the weighted averages.
    # TODO: need a more refined way of computing weights when block= is set.
    ncells <- table(groups)
    if (is.null(weight.fun)) {
        weight.fun <- function(x) log10(pmin(x, 100))
    }
    w <- weight.fun(ncells)
    names(w) <- names(ncells)

    summary.effects <- vector("list", length(output.effect[[1]]))
    for (i in seq_along(output.effect[[1]])) {
        current.out <- list()

        for (out in names(output.effect)) {
            effect.df <- output.effect[[out]][[i]]
            effect.mat <- as.matrix(effect.df)

            current.out[[paste0("mean.", out)]] <- rowWeightedMeans(effect.mat, w=w[colnames(effect.mat)], na.rm=TRUE)
            current.out[[paste0("min.", out)]] <- rowMins(effect.mat, na.rm=TRUE)
            current.out[[paste0("med.", out)]] <- rowWeightedMedians(effect.mat, w=w[colnames(effect.mat)], na.rm=TRUE)
            current.out[[paste0("max.", out)]] <- rowMaxs(effect.mat, na.rm=TRUE)
        
            if (full.stats) {
                current.out[[paste0("full.", out)]] <- I(effect.df)
            }
        }

        df <- do.call(DataFrame, current.out)
        rownames(df) <- rownames(output.effect[[out]][[1]])

        summary.effects[[i]] <- df
    }

    names(summary.effects) <- names(output.effect[[1]])

    # Throwing in the usual stats.
    row.data <- summaryMarkerStats(x, groups, BPPARAM=BPPARAM)
    summary.effects <- .add_row_data(summary.effects, row.data, match.names=FALSE)

    List(summary.effects)
}


