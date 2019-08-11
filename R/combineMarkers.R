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
#' @param pval.type A string specifying the type of combined p-value to be computed, i.e., Simes' (\code{"any"}) or IUT (\code{"all"}).
#' @param log.p.in A logical scalar indicating if the p-values in \code{de.lists} were log-transformed.
#' @param log.p.out A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param output.field A string specifying the prefix of the field names containing the effect sizes.
#' Defaults to \code{"stats"} if \code{full.stats=TRUE}, otherwise it is set to \code{effect.field}.
#' @param full.stats A logical scalar indicating whether all statistics in \code{de.lists} should be stored in the output for each pairwise comparison.
#' @param sorted Logical scalar indicating whether each output DataFrame should be sorted by a statistic relevant to \code{pval.type}.
#' 
#' @return
#' A named \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains the consolidated results for the cluster of the same name.
#' 
#' Within each DataFrame (say, the DataFrame for cluster X), rows correspond to genes with the fields:
#' \describe{
#' \item{\code{Top}:}{Integer, the minimum rank across all pairwise comparisons.
#' Only reported if \code{pval.type="any"}.}
#' \item{\code{p.value}:}{Numeric, the p-value across all comparisons if \code{log.p.out=FALSE}.
#' This is a Simes' p-value if \code{pval.type="any"}, otherwise it is an IUT p-value.}
#' \item{\code{log.p.value}:}{Numeric, the log-transformed version of \code{p.value} if \code{log.p.out=TRUE}.
#' Natural logarithms are reported.}
#' \item{\code{FDR}:}{Numeric, the BH-adjusted p-value for each gene if \code{log.p.out=FALSE}.}
#' \item{\code{log.FDR}:}{Numeric, the log-transformed adjusted p-value for each gene if \code{log.p.out=TRUE}.
#' Natural logarithms are reported.}
#' \item{\code{logFC.Y}:}{Numeric field present for every other cluster Y in \code{clusters}.
#' It contains the effect size of the comparison of X to Y:
#' \itemize{
#' \item Log-fold changes, usually in base 2, from \code{\link{pairwiseTTests}}.
#' \item Area under the curve from \code{\link{pairwiseWilcox}}.
#' \item Log2-fold changes in the expressing proportion from \code{\link{pairwiseBinom}}.
#' }
#' The exact name of this field depends on \code{output.field}, see Examples.
#' It is only reported when \code{full.stats=FALSE}.}
#' \item{\code{stats.Y}:}{DataFrame for every other cluster Y in \code{clusters}, returned when \code{full.stats=TRUE}.
#' This contains the same fields in the corresponding entry of \code{de.lists} for the X versus Y comparison.
#' The name of this field can be altered by setting \code{output.field}.}
#' }
#' If \code{sorted=TRUE}, genes are ranked by the \code{Top} column (if available) and then the \code{p.value} column.
#' 
#' The DataFrames themselves are sorted according to the order of cluster IDs in \code{pairs[,1]}.
#' The \code{logFC.Y} columns are sorted according to the order of cluster IDs in \code{pairs[,2]} within the corresponding level of the first cluster.
#' 
#' Note that DataFrames are only created for clusters present in \code{pairs[,1]}.
#' Clusters unique to \code{pairs[,2]} will only be present within each DataFrame as Y.
#' 
#' @details
#' An obvious strategy to characterizing differences between clusters is to look for genes that are differentially expressed (DE) between them.
#' However, this entails a number of comparisons between all pairs of clusters to comprehensively identify genes that define each cluster.
#' For all pairwise comparisons involving a single cluster, we would like to consolidate the DE results into a single list of candidate marker genes.
#' This is the intention of the \code{combineMarkers} function.
#' DE statistics from any testing regime can be supplied to this function - see the Examples for how this is done with t-tests from \code{\link{pairwiseTTests}}.
#' 
#' @section Consolidating p-values into a ranking:
#' By default, each table is sorted by the \code{Top} value when \code{pval.type="any"}.
#' This is the minimum rank across all pairwise comparisons for each gene, and specifies the size of the candidate marker set.
#' Taking all rows with \code{Top} values less than or equal to X will yield a marker set containing the top X genes (ranked by significance) from each pairwise comparison.
#' The marker set for each cluster allows it to be distinguished from every other cluster based on the differential expression of at least one gene.
#' 
#' To demonstrate, let us define a marker set with an X of 1 for a given cluster.
#' The set of genes with \code{Top <= 1} will contain the top gene from each pairwise comparison to every other cluster.
#' If X is instead, say, 5, the set will consist of the \emph{union} of the top 5 genes from each pairwise comparison.
#' Obviously, multiple genes can have the same \code{Top} as different genes may have the same rank across different pairwise comparisons.
#' Conversely, the marker set may be smaller than the product of \code{Top} and the number of other clusters, as the same gene may be shared across different comparisons.
#' 
#' This approach does not explicitly favour genes that are uniquely expressed in a cluster.
#' Such a strategy is often too stringent, especially in cases involving overclustering or cell types defined by combinatorial gene expression.
#' However, if \code{pval.type="all"}, the null hypothesis is that the gene is not DE in all contrasts, and the IUT p-value is computed for each gene.
#' This yields a \code{IUT.p} field instead of a \code{Top} field in the output table.
#' Ranking based on the IUT p-value will focus on genes that are uniquely DE in that cluster.
#' 
#' @section Correcting for multiple testing:
#' When \code{pval.type="any"}, a combined p-value is calculated by consolidating p-values across contrasts for each gene using Simes' method.
#' This represents the evidence against the null hypothesis is that the gene is not DE in any of the contrasts.
#' The BH method is then applied on the combined p-values across all genes to obtain the \code{FDR} field.
#' The same procedure is done with \code{pval.type="all"}, but using the IUT p-values across genes instead.
#' 
#' If \code{log.p=TRUE}, log-transformed p-values and FDRs will be reported.
#' This may be useful in over-powered studies with many cells, where directly reporting the raw p-values would result in many zeroes due to the limits of machine precision.
#' 
#' Note that the reported FDRs are intended only as a rough measure of significance.
#' Properly correcting for multiple testing is not generally possible when \code{clusters} is determined from the same \code{x} used for DE testing.
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
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(sce)), centers=3)
#' clusters <- paste0("Cluster", kout$cluster)
#' 
#' out <- pairwiseTTests(logcounts(sce), clusters=clusters)
#' comb <- combineMarkers(out$statistics, out$pairs)
#' comb[["Cluster1"]]
#' 
#' out <- pairwiseWilcox(logcounts(sce), clusters=clusters)
#' comb <- combineMarkers(out$statistics, out$pairs, effect.field="AUC")
#' comb[["Cluster2"]]
#' 
#' out <- pairwiseBinom(logcounts(sce), clusters=clusters)
#' comb <- combineMarkers(out$statistics, out$pairs)
#' comb[["Cluster3"]]
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame List
#' @importFrom stats p.adjust
#' @importFrom BiocGenerics cbind
#' @importFrom methods as
combineMarkers <- function(de.lists, pairs, pval.field="p.value", effect.field="logFC", 
    pval.type=c("any", "all"), log.p.in=FALSE, log.p.out=log.p.in, 
    output.field=NULL, full.stats=FALSE, sorted=TRUE)
{
    if (length(de.lists)!=nrow(pairs)) {
        stop("'nrow(pairs)' must be equal to 'length(de.lists)'")
    }
    if (is.null(output.field)) {
        output.field <- if (full.stats) "stats" else effect.field
    }

    pval.type <- match.arg(pval.type)
    method <- switch(pval.type, any="simes", all="berger")

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
    output <- vector("list", length(by.first))
    names(output) <- names(by.first)

    for (host in names(by.first)) {
        chosen <- by.first[[host]]
        targets <- pairs[chosen, 2]
        cur.stats <- de.lists[chosen]

        keep <- !is.na(targets)
        targets <- targets[keep]
        cur.stats <- cur.stats[keep]

        all.p <- lapply(cur.stats, "[[", i=pval.field)
        pval <- do.call(combinePValues, c(all.p, list(method=method, log.p=log.p.in)))
        preamble <- DataFrame(row.names=gene.names)

        # Determining rank.
        if (pval.type=="any") {
            rank.out <- .rank_top_genes(all.p)
            min.rank <- rank.out$rank
            min.p <- rank.out$value
            gene.order <- order(min.rank, min.p)
            preamble$Top <- min.rank
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
        preamble[[paste0(prefix, "p.value")]] <- pval 
        preamble[[paste0(prefix, "FDR")]] <- corrected 

        # Saving effect sizes or all statistics.
        if (full.stats) {
            cur.stats <- lapply(cur.stats, FUN=function(x) { I(as(x, Class="DataFrame")) })
            stat.df <- do.call(DataFrame, c(cur.stats, list(check.names=FALSE)))
        } else {
            all.effects <- lapply(cur.stats, "[[", i=effect.field)
            stat.df <- DataFrame(all.effects)
        }
        colnames(stat.df) <- sprintf("%s.%s", output.field, targets)

        # Producing the output object.
        marker.set <- cbind(preamble, stat.df)
        if (sorted) {
            marker.set <- marker.set[gene.order,,drop=FALSE]
        }
        output[[host]] <- marker.set
    }

    as(output, "List")
}

.rank_top_genes <- function(metrics) 
# This computes the rank and the minimum metric for each gene.
{
    ncon <- length(metrics)
    ngenes <- if (ncon) length(metrics[[1]]) else 0L
    min.rank <- min.val <- rep(NA_integer_, ngenes)

    for (con in seq_len(ncon)) { 
        cur.val <- metrics[[con]]
        cur.rank <- rank(cur.val, ties.method="first", na.last="keep")
        min.rank <- pmin(min.rank, cur.rank, na.rm=TRUE)
        min.val <- pmin(min.val, cur.val, na.rm=TRUE)
    }
    
    list(rank=min.rank, value=min.val)
}

