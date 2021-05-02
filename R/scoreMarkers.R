#' Score markers
#'
#' Compute various summary scores for potential marker genes to distinguish between groups of cells.
#'
#' @param x A matrix-like object containing log-expression values, with genes in rows and cells in columns.
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
#' Each DataFrame corresponds to a group and each row corresponds to a marker in \code{x}.
#' See Details for information about the individual columns.
#'
#' @details
#' Compared to \code{\link{findMarkers}}, this function represents a simpler and more intuitive summary of the differences between the groups.
#' We do this by realizing that the p-values for these types of comparisons are largely meaningless;
#' individual cells are not meaningful units of experimental replication, while the groups themselves are defined from the data.
#' Thus, by discarding the p-values, we can simplify our marker selection by focusing only on the effect sizes between groups.
#' 
#' Here, the strategy is to perform pairwise comparisons between each pair of groups to obtain various effect sizes.
#' We consider three different effects:
#' \itemize{
#' \item \code{logFC.cohen}, the standardized log-fold change.
#' This is the difference in the mean log-expression for each group scaled by the root of the pooled variance across the groups.
#' The standardization adjusts for differences in variances between comparisons and is analogous to the calculation of the t-statistic.
#' \item \code{AUC}, the area under the curve.
#' This is the probability that a randomly chosen observation in one group is greater than a randomly chosen observation in the other group.
#' The AUC is closely related to the U-statistic used in the Wilcoxon rank sum test.
#' \item \code{logFC.detected}, the log-fold change in the proportion of cells with detected expression between groups.
#' }
#' 
#' For each group \eqn{X}, we then consider all pairwise comparisons between \eqn{X} and other groups. 
#' The statistics from all of these pairwise comparisons are summarized into a few metrics in the output DataFrame for \eqn{X}.
#' For example, for the AUC, we have:
#' \itemize{
#' \item \code{mean.AUC}, the mean AUC across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to the average of the other groups.
#' \item \code{median.AUC}, the median AUC across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to most other groups.
#' \item \code{min.AUC}, the minimum AUC across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to all other groups.
#' \item \code{max.AUC}, the minimum AUC across all pairwise comparisons involving \eqn{X}.
#' A positive value indicates that the gene is upregulated in \eqn{X} compared to at least one other group.
#' }
#' Similar columns are generated for the other effects.
#'
#' In practice, the mean and the median metrics are obtained by weighting each comparison according to the number of cells involved.
#' This ensures that statistics from comparisons with very few cells do not skew the summaries.
#' If \code{weight.fun=NULL}, the default weight is defined as the log10-transformed number of cells;
#' this favors comparisons with more cells but ensures that comparisons involving very large groups do not dominate the summary.
#' 
#' If \code{full.stats=TRUE}, an extra \code{full.*} column is returned in the DataFrame.
#' This contains a nested DataFrame with number of columns equal to the number of other groups.
#' Each column contains the statistic from the comparison between \eqn{X} and the other group.
#' 
#' In addition, we report the mean log-expression of \eqn{X} as well as the grand mean of mean log-expression values across all other groups.
#' The same is done for the proportion of cells with detectable expression.
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

    for (tt in c("wilcox", "binom", "t")) {
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
        weight.fun <- log10
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


