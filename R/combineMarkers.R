#' Combine pairwise DE results into a marker list
#' 
#' Combine multiple pairwise differential expression comparisons between groups or clusters into a single ranked list of markers for each cluster.
#' 
#' @param de.lists A list-like object where each element is a data.frame or \linkS4class{DataFrame}.
#' Each element should represent the results of a pairwise comparison between two groups/clusters,
#' in which each row should contain the statistics for a single gene/feature.
#' Rows should be named by the feature name in the same order for all elements.
#' @param pairs A matrix, data.frame or \linkS4class{DataFrame} with two columns and number of rows equal to the length of \code{de.lists}.
#' Each row should specify the pair of clusters being compared for the corresponding element of \code{de.lists}.
#' @param pval.field A string specifying the column name of each element of \code{de.lists} that contains the p-value.
#' @param effect.field A string specifying the column name of each element of \code{de.lists} that contains the effect size.
#' If \code{NULL}, effect sizes are not reported in the output.
#' @param pval.type A string specifying how p-values are to be combined across pairwise comparisons for a given group/cluster.
#' @param min.prop Numeric scalar specifying the minimum proportion of significant comparisons per gene,
#' Defaults to 0.5 when \code{pval.type="some"}, otherwise defaults to zero.
#' @param log.p.in A logical scalar indicating if the p-values in \code{de.lists} were log-transformed.
#' @param log.p.out A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param output.field A string specifying the prefix of the field names containing the effect sizes.
#' Defaults to \code{"stats"} if \code{full.stats=TRUE}, otherwise it is set to \code{effect.field}.
#' @param full.stats A logical scalar indicating whether all statistics in \code{de.lists} should be stored in the output for each pairwise comparison.
#' @param sorted Logical scalar indicating whether each output DataFrame should be sorted by a statistic relevant to \code{pval.type}.
#' @param flatten Logical scalar indicating whether the individual effect sizes should be flattened in the output DataFrame.
#' If \code{FALSE}, effect sizes are reported as a nested matrix for easier programmatic use.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes.
#' 
#' @return
#' A named \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains the consolidated marker statistics for each gene (row) for the cluster of the same name.
#' The DataFrame for cluster \eqn{X} contains the fields:
#' \describe{
#' \item{\code{Top}:}{Integer, the minimum rank across all pairwise comparisons - see \code{?\link{computeMinRank}} for details.
#' This is only reported if \code{pval.type="any"}.}
#' \item{\code{p.value}:}{Numeric, the combined p-value across all comparisons if \code{log.p.out=FALSE}.}
#' \item{\code{FDR}:}{Numeric, the BH-adjusted p-value for each gene if \code{log.p.out=FALSE}.}
#' \item{\code{log.p.value}:}{Numeric, the (natural) log-transformed version p-value.
#' Replaces the \code{p.value} field if \code{log.p.out=TRUE}.}
#' \item{\code{log.FDR}:}{Numeric, the (natural) log-transformed adjusted p-value.
#' Replaces the \code{FDR} field if \code{log.p.out=TRUE}.}
#' \item{\code{summary.<OUTPUT>}:}{Numeric, named by replacing \code{<OUTPUT>} with \code{output.field}.
#' This contains the summary effect size, obtained by combining effect sizes from all pairwise comparison into a single value.
#' Only reported when \code{effect.field} is not \code{NULL}.}
#' \item{\code{<OUTPUT>.Y}:}{Comparison-specific statistics, named by replacing \code{<OUTPUT>} with \code{output.field}.
#' One of these fields is present for every other cluster \eqn{Y} in \code{clusters} and contains statistics for the comparison of \eqn{X} to \eqn{Y}.
#' If \code{full.stats=FALSE}, each field is numeric and contains the effect size of the comparison of \eqn{X} over \eqn{Y}.
#' Otherwise, each field is a nested DataFrame containing the full statistics for that comparison (i.e., the same asthe corresponding entry of \code{de.lists}).
#' Only reported if \code{flatten=FALSE} and (for \code{full.stats=FALSE}) if \code{effect.field} is not \code{NULL}.
#' }
#' \item{\code{each.<OUTPUT>}:}{A nested DataFrame of comparison-specific statistics, named by replacing \code{<OUTPUT>} with \code{output.field}.
#' If \code{full.stats=FALSE}, one column is present for every other cluster \eqn{Y} in \code{clusters} and contains the effect size of the comparison of \eqn{X} to \eqn{Y}.
#' Otherwise, each column contains another nested DataFrame containing the full set of statistics for that comparison.
#' Only reported if \code{flatten=FALSE} and (for \code{full.stats=FALSE}) if \code{effect.field} is not \code{NULL}.
#' }
#' }
#' 
#' @details
#' An obvious strategy to characterizing differences between clusters is to look for genes that are differentially expressed (DE) between them.
#' However, this entails a number of comparisons between all pairs of clusters to comprehensively identify genes that define each cluster.
#' For all pairwise comparisons involving a single cluster, we would like to consolidate the DE results into a single list of candidate marker genes.
#' Doing so is the purpose of the \code{combineMarkers} function.
#'
#' DE statistics from any testing regime can be supplied to this function - see the Examples for how this is done with t-tests from \code{\link{pairwiseTTests}}.
#' The effect size field in the output will vary according to the type of input statistics, for example:
#' \itemize{
#' \item \code{logFC.Y} from \code{\link{pairwiseTTests}}, containing log-fold changes in mean expression (usually in base 2).
#' \item \code{AUC.Y} from \code{\link{pairwiseWilcox}}, containing the area under the curve, i.e., the concordance probability. 
#' \item \code{logFC.Y} from \code{\link{pairwiseBinom}}, containing log2-fold changes in the expressing proportion.
#' }
#' 
#' @section Consolidating with DE against any other cluster:
#' When \code{pval.type="any"}, each DataFrame is sorted by the min-rank in the \code{Top} column.
#' Taking all rows with \code{Top} values less than or equal to \eqn{T} yields a marker set containing the top \eqn{T} genes (ranked by significance) from each pairwise comparison.
#' This guarantees the inclusion of genes that can distinguish between any two clusters.
#' Also see \code{?\link{computeMinRank}} for more details on the rationale behind this metric.
#' 
#' For each gene and cluster, the summary effect size is defined as the effect size from the pairwise comparison with the lowest p-value.
#' The combined p-value is computed by applying Simes' method to all p-values.
#' Neither of these values are directly used for ranking and are only reported for the sake of the user.
#'
#' @section Consolidating with DE against all other clusters:
#' If \code{pval.type="all"}, the null hypothesis is that the gene is not DE in all contrasts.
#' A combined p-value for each gene is computed using Berger's intersection union test (IUT).
#' Ranking based on the IUT p-value will focus on genes that are DE in that cluster compared to \emph{all} other clusters.
#' This strategy is particularly effective when dealing with distinct clusters that have a unique expression profile.
#' In such cases, it yields a highly focused marker set that concisely captures the differences between clusters.
#' 
#' However, it can be too stringent if the cluster's separation is driven by combinations of gene expression.
#' For example, consider a situation involving four clusters expressing each combination of two marker genes A and B.
#' With \code{pval.type="all"}, neither A nor B would be detected as markers as it is not uniquely defined in any one cluster.
#' This is especially detrimental with overclustering where an otherwise acceptable marker is discarded if it is not DE between two adjacent clusters.
#'
#' For each gene and cluster, the summary effect size is defined as the effect size from the pairwise comparison with the \emph{largest} p-value.
#' This reflects the fact that, with this approach, a gene is only as significant as its weakest DE.
#' Again, this value is not directly used for ranking and are only reported for the sake of the user.
#'
#' @section Consolidating with DE against some other clusters:
#' The \code{pval.type="some"} setting serves as a compromise between \code{"all"} and \code{"any"}.
#' A combined p-value is calculated by taking the middlemost value of the Holm-corrected p-values for each gene.
#' (By default, this the median for odd numbers of contrasts and one-after-the-median for even numbers, but the exact proportion can be changed by setting \code{min.prop} - see \code{?\link{combineParallelPValues}}.)
#' Here, the null hypothesis is that the gene is not DE in at least half of the contrasts.
#' 
#' Genes are then ranked by the combined p-value.
#' The aim is to provide a more focused marker set without being overly stringent, though obviously it loses the theoretical guarantees of the more extreme settings.
#' For example, there is no guarantee that the top set contains genes that can distinguish a cluster from any other cluster, which would have been possible with \code{pval.type="any"}.
#'
#' For each gene and cluster, the summary effect size is defined as the effect size from the pairwise comparison with the \code{min.prop}-smallest p-value.
#' This mirrors the p-value calculation but, again, is reported only for the benefit of the user.
#'
#' @section Consolidating against some other clusters, rank-style:
#' A slightly different flavor of the \dQuote{some cluster} approach is achieved by setting \code{method="any"} with \code{min.prop} set to some positive value in (0, 1).
#' A gene will only be high-ranked if it is among the top-ranked genes in at least \code{min.prop} of the pairwise comparisons.
#' For example, if \code{min.prop=0.3}, any gene with a value of \code{Top} less than or equal to 5 will be in the top 5 DEGs of at least 30\% of the comparisons.
#' 
#' This method increases the stringency of the \code{"any"} setting in a safer manner than \code{pval.type="some"}.
#' Specifically, we avoid comparing p-values across pairwise comparisons, which can be problematic if there are power differences across comparisons, e.g., due to differences in the number of cells across the other clusters.
#'
#' Note that the value of \code{min.prop} does not affect the combined p-value and summary effect size calculations for \code{pval.type="any"}.
#'
#' @section Correcting for multiple testing:
#' The BH method is then applied on the consolidated p-values across all genes to obtain the \code{FDR} field.
#' The reported FDRs are intended only as a rough measure of significance.
#' Properly correcting for multiple testing is not generally possible when \code{clusters} is determined from the same \code{x} used for DE testing.
#' 
#' If \code{log.p=TRUE}, log-transformed p-values and FDRs will be reported.
#' This may be useful in over-powered studies with many cells, where directly reporting the raw p-values would result in many zeroes due to the limits of machine precision.
#' 
#' @section Ordering of the output:
#' \itemize{
#' \item Within each DataFrame, if \code{sorted=TRUE}, genes are ranked by the \code{Top} column if available and the \code{p.value} (or \code{log.p.value}) if not.
#' Otherwise, the input order of the genes is preserved.
#' \item For the DataFrame corresponding to cluster \eqn{X}, the \code{<OUTPUT>.Y} columns are sorted according to the order of cluster IDs in \code{pairs[,2]} for all rows where \code{pairs[,1]} is \eqn{X}.
#' \item In the output List, the DataFrames themselves are sorted according to the order of cluster IDs in \code{pairs[,1]}.
#' Note that DataFrames are only created for clusters present in \code{pairs[,1]}.
#' Clusters unique to \code{pairs[,2]} will only be present within a DataFrame as \eqn{Y}.
#' }
#' 
#' @seealso
#' \code{\link{pairwiseTTests}} and \code{\link{pairwiseWilcox}}, for functions that can generate \code{de.lists} and \code{pairs}.
#' 
#' \code{\link{findMarkers}}, which automatically performs \code{combineMarkers} on the t-test or Wilcoxon test results.
#' 
#' @references
#' Simes RJ (1986). 
#' An improved Bonferroni procedure for multiple tests of significance. 
#' \emph{Biometrika} 73:751-754.
#' 
#' Berger RL and Hsu JC (1996). 
#' Bioequivalence trials, intersection-union tests and equivalence confidence sets.
#' \emph{Statist. Sci.} 11, 283-319.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(sce)), centers=3)
#' clusters <- paste0("Cluster", kout$cluster)
#' 
#' out <- pairwiseTTests(logcounts(sce), groups=clusters)
#' comb <- combineMarkers(out$statistics, out$pairs)
#' comb[["Cluster1"]]
#' 
#' out <- pairwiseWilcox(logcounts(sce), groups=clusters)
#' comb <- combineMarkers(out$statistics, out$pairs, effect.field="AUC")
#' comb[["Cluster2"]]
#' 
#' out <- pairwiseBinom(logcounts(sce), groups=clusters)
#' comb <- combineMarkers(out$statistics, out$pairs)
#' comb[["Cluster3"]]
#'
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom metapod combineParallelPValues
combineMarkers <- function(de.lists, pairs, pval.field="p.value", effect.field="logFC", 
    pval.type=c("any", "some", "all"), min.prop=NULL, log.p.in=FALSE, log.p.out=log.p.in, 
    output.field=NULL, full.stats=FALSE, sorted=TRUE, flatten=TRUE, BPPARAM=SerialParam())
{
    if (length(de.lists)!=nrow(pairs)) {
        stop("'nrow(pairs)' must be equal to 'length(de.lists)'")
    }
    if (is.null(output.field)) {
        output.field <- if (full.stats) "stats" else effect.field
    }
    report.effects <- !is.null(effect.field)

    pval.type <- match.arg(pval.type)
    method <- switch(pval.type, any="simes", some="holm-min", all="berger")
    if (is.null(min.prop))  {
        min.prop <- if (pval.type=="any") 0 else 0.5
    }

    # Checking that all genes are the same across lists.
    gene.names <- NULL
    for (x in seq_along(de.lists)) {
        current <- de.lists[[x]]
        curnames <- rownames(current)

        if (is.null(gene.names)) {
            gene.names <- curnames
        } else if (!identical(gene.names, curnames)) {
            stop("row names should be the same for all elements of 'de.lists'")
        }
    }

    # Processing by the first element of each pair.
    first.fac <- factor(pairs[,1], unique(pairs[,1]))
    by.first <- split(seq_along(de.lists), first.fac, drop=TRUE)

    output <- bplapply(by.first, FUN=.combine_markers_internal,
        pairs=pairs, de.lists=de.lists, method=method, gene.names=gene.names, report.effects=report.effects,
        pval.type=pval.type, effect.field=effect.field, min.prop=min.prop, log.p.in=log.p.in, log.p.out=log.p.out, 
        pval.field=pval.field, output.field=output.field, full.stats=full.stats, sorted=sorted, flatten=flatten, 
        BPPARAM=BPPARAM)

    SimpleList(output)
}

#' @importFrom methods as
#' @importFrom stats p.adjust
#' @importFrom BiocGenerics cbind
#' @importFrom S4Vectors DataFrame I
.combine_markers_internal <- function(chosen, pairs, de.lists, method, gene.names, report.effects, 
    pval.type, effect.field, min.prop, log.p.in, log.p.out, pval.field, output.field, full.stats,
    sorted, flatten, BPPARAM=SerialParam())
{
    targets <- pairs[chosen, 2]
    cur.stats <- de.lists[chosen]

    keep <- !is.na(targets)
    targets <- targets[keep]
    cur.stats <- cur.stats[keep]

    all.p <- lapply(cur.stats, "[[", i=pval.field)
    pval <- combineParallelPValues(all.p, method=method, log.p=log.p.in, min.prop=min.prop)$p.value
    marker.set <- DataFrame(row.names=gene.names)

    # Determining rank.
    if (pval.type=="any") {
        ranked <- lapply(all.p, rank, ties.method="first", na.last="keep")
        min.rank <- compute_Top_statistic_from_ranks(ranked, min.prop)
        gene.order <- order(min.rank) 
        marker.set$Top <- min.rank
    } else {
        gene.order <- order(pval)
    }

    # Correcting for multiple testing. We try to preserve the log-ness as long as we can,
    # to avoid underflow upon exp()'ing that could be avoided by correction.
    if (log.p.in) {
        corrected <- .logBH(pval)
    } else {
        corrected <- p.adjust(pval, method="BH")
    }
    if (log.p.out!=log.p.in) {
        transFUN <- if (log.p.out) log else exp
        pval <- transFUN(pval)
        corrected <- transFUN(corrected)
    }
    
    prefix <- if (log.p.out) "log." else ""
    marker.set[[paste0(prefix, "p.value")]] <- pval 
    marker.set[[paste0(prefix, "FDR")]] <- corrected 

    if (report.effects) {
        all.effects <- lapply(cur.stats, "[[", i=effect.field)
        marker.set[[paste0("summary.", output.field)]] <- .choose_effect_size(all.p, all.effects, pval.type, min.prop)
    }

    # Saving effect sizes or all statistics.
    if (full.stats || report.effects) {
        if (full.stats) {
            cur.stats <- lapply(cur.stats, FUN=function(x) { I(DataFrame(x)) })
        } else {
            cur.stats <- all.effects
        }
        stat.df <- do.call(DataFrame, c(cur.stats, list(row.names=gene.names)))

        if (flatten) {
            colnames(stat.df) <- sprintf("%s.%s", output.field, targets)
            marker.set <- cbind(marker.set, stat.df)
        } else {
            colnames(stat.df) <- as.character(targets)
            marker.set[[paste0("each.", output.field)]] <- stat.df
        }
    }

    if (sorted) {
        marker.set <- marker.set[gene.order,,drop=FALSE]
    }
    marker.set
}

.choose_effect_size <- function(all.p, all.effects, pval.type, min.prop) {
    if (!length(all.p)) {
        return(numeric(0))
    }

    ngenes <- length(all.p[[1]])
    chosen <- rep(NA_real_, ngenes)

    if (pval.type=="all") {
        # Don't start from worst <- all.p[[i]] as this would
        # require another layer of protection against NA's.
        worst <- rep(-Inf, ngenes)
        for (i in seq_along(all.p)) {
            worse <- all.p[[i]] > worst & !is.na(all.p[[i]])
            chosen[worse] <- all.effects[[i]][worse]
            worst[worse] <- all.p[[i]][worse]
        }

    } else if (pval.type=="any") {
        best <- rep(Inf, ngenes)
        for (i in seq_along(all.p)) {
            better <- all.p[[i]] < best & !is.na(all.p[[i]])
            chosen[better] <- all.effects[[i]][better]
            best[better] <- all.p[[i]][better]
        }

    } else {
        chosen <- choose_middle_effect_size(all.p, all.effects, min.prop)
    }

    chosen
}
