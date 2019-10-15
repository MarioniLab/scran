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
#' @param log.p.in A logical scalar indicating if the p-values in \code{de.lists} were log-transformed.
#' @param log.p.out A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param output.field A string specifying the prefix of the field names containing the effect sizes.
#' Defaults to \code{"stats"} if \code{full.stats=TRUE}, otherwise it is set to \code{effect.field}.
#' @param full.stats A logical scalar indicating whether all statistics in \code{de.lists} should be stored in the output for each pairwise comparison.
#' @param sorted Logical scalar indicating whether each output DataFrame should be sorted by a statistic relevant to \code{pval.type}.
#' 
#' @return
#' A named \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains the consolidated marker statistics for each gene (row) for the cluster of the same name.
#' The DataFrame for cluster \eqn{X} contains the fields:
#' \describe{
#' \item{\code{Top}:}{Integer, the minimum rank across all pairwise comparisons.
#' This is only reported if \code{pval.type="any"}.}
#' \item{\code{p.value}:}{Numeric, the p-value across all comparisons if \code{log.p.out=FALSE}.
#' This is a Simes' p-value if \code{pval.type="any"} or \code{"some"}, otherwise it is an IUT p-value.}
#' \item{\code{FDR}:}{Numeric, the BH-adjusted p-value for each gene if \code{log.p.out=FALSE}.}
#' \item{\code{log.p.value}:}{Numeric, the (natural) log-transformed version p-value.
#' Replaces the \code{p.value} field if \code{log.p.out=TRUE}.}
#' \item{\code{log.FDR}:}{Numeric, the (natural) log-transformed adjusted p-value.
#' Replaces the \code{FDR} field if \code{log.p.out=TRUE}.}
#' \item{\code{<OUTPUT>.Y}:}{Numeric, where the value of \code{<OUTPUT>} is set to \code{output.field}.
#' One of these fields is present for every other cluster \eqn{Y} in \code{clusters}.
#' It contains the effect size of the comparison of \eqn{X} to \eqn{Y}, with typical values depending on the chosen method:
#' \itemize{
#' \item \code{logFC.Y} from \code{\link{pairwiseTTests}}, containing log-fold changes in mean expression (usually in base 2).
#' \item \code{AUC.Y} from \code{\link{pairwiseWilcox}}, containing the area under the curve, i.e., the concordance probability. 
#' \item \code{logFC.Y} from \code{\link{pairwiseBinom}}, containing log2-fold changes in the expressing proportion.
#' }
#' Only reported when \code{full.stats=FALSE} and \code{effect.field} is not \code{NULL}.
#' }
#' \item{\code{stats.y}:}{A nested DataFrame containing the same fields in the corresponding entry of \code{de.lists} for the \eqn{X} versus \eqn{Y} comparison.
#' Only reported when \code{full.stats=TRUE}.}
#' }
#' 
#' Ordering rules are:
#' \itemize{
#' \item Within each DataFrame, if \code{sorted=TRUE}, genes are ranked by the \code{Top} column (if available) and then the \code{p.value} column.
#' Otherwise, the input order of the genes is preserved.
#' \item For the DataFrame corresponding to cluster X, the \code{<OUTPUT>.Y} columns are sorted according to the order of cluster IDs in \code{pairs[,2]} for all rows where \code{pairs[,1]} is X.
#' \item In the output List, the DataFrames themselves are sorted according to the order of cluster IDs in \code{pairs[,1]}.
#' Note that DataFrames are only created for clusters present in \code{pairs[,1]}.
#' Clusters unique to \code{pairs[,2]} will only be present within a DataFrame as Y.
#' }
#' 
#' @details
#' An obvious strategy to characterizing differences between clusters is to look for genes that are differentially expressed (DE) between them.
#' However, this entails a number of comparisons between all pairs of clusters to comprehensively identify genes that define each cluster.
#' For all pairwise comparisons involving a single cluster, we would like to consolidate the DE results into a single list of candidate marker genes.
#' Doing so is the purpose of the \code{combineMarkers} function.
#' DE statistics from any testing regime can be supplied to this function - see the Examples for how this is done with t-tests from \code{\link{pairwiseTTests}}.
#' 
#' @section Consolidating with DE against any other cluster:
#' By default, each DataFrame is sorted by the \code{Top} value when \code{pval.type="any"}.
#' (For genes with the same \code{Top}, ranking is performed based on the Simes combined p-value - see the comments for \code{pval.type="some"}.)
#' Taking all rows with \code{Top} values less than or equal to T yields a marker set containing the top T genes (ranked by significance) from each pairwise comparison.
#' This guarantees the inclusion of genes that can distinguish between any two clusters.
#' 
#' To demonstrate, let us define a marker set with an T of 1 for a given cluster.
#' The set of genes with \code{Top <= 1} will contain the top gene from each pairwise comparison to every other cluster.
#' If T is instead, say, 5, the set will consist of the \emph{union} of the top 5 genes from each pairwise comparison.
#' Obviously, multiple genes can have the same \code{Top} as different genes may have the same rank across different pairwise comparisons.
#' Conversely, the marker set may be smaller than the product of \code{Top} and the number of other clusters, as the same gene may be shared across different comparisons..
#' 
#' This approach does not explicitly favour genes that are uniquely expressed in a cluster.
#' Rather, it focuses on combinations of genes that - together - drive separation of a cluster from the others.
#' This is more general and robust but tends to yield a less focused marker set compared to the other \code{pval.type} settings.
#'
#' @section Consolidating with DE against all other clusters:
#' If \code{pval.type="all"}, the null hypothesis is that the gene is not DE in all contrasts, and the IUT p-value is computed for each gene.
#' Ranking based on the IUT p-value will focus on genes that are DE in that cluster compared to \emph{all} other clusters.
#' This strategy is particularly effective when dealing with distinct clusters that have a unique expression profile.
#' In such cases, it yields a highly focused marker set that concisely captures the differences between clusters.
#' 
#' However, it can be too stringent if the cluster's separation is driven by combinations of gene expression.
#' For example, consider a situation involving four clusters expressing each combination of two marker genes A and B.
#' With \code{pval.type="all"}, neither A nor B would be detected as markers as it is not uniquely defined in any one cluster.
#' This is especially detrimental with overclustering where an otherwise acceptable marker is discarded if it is not DE between two adjacent clusters.
#'
#' @section Consolidating with DE against some other clusters:
#' The \code{pval.type="some"} setting serves as a compromise between \code{"all"} and \code{"any"}.
#' A combined p-value is calculated by taking the middlemost value of the Holm-corrected p-values for each gene.
#' (This is the median for odd numbers of contrasts and one-after-the-median for even numbers, see \code{?\link{combinePValues}}.)
#' Here, the null hypothesis is that the gene is not DE in at least half of the contrasts.
#' 
#' Genes are then ranked by the combined p-value.
#' The aim is to provide a more focused marker set without being overly stringent, though obviously it loses the theoretical guarantees of the more extreme settings.
#' For example, there is no guarantee that the top set contains genes that can distinguish a cluster from any other cluster, which would have been possible with \code{"any"}.
#'
#' @section Correcting for multiple testing:
#' The BH method is then applied on the Simes/IUT p-values across all genes to obtain the \code{FDR} field.
#' The reported FDRs are intended only as a rough measure of significance.
#' Properly correcting for multiple testing is not generally possible when \code{clusters} is determined from the same \code{x} used for DE testing.
#' 
#' If \code{log.p=TRUE}, log-transformed p-values and FDRs will be reported.
#' This may be useful in over-powered studies with many cells, where directly reporting the raw p-values would result in many zeroes due to the limits of machine precision.
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
#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame 
#' @importFrom stats p.adjust
#' @importFrom BiocGenerics cbind
#' @importFrom methods as
combineMarkers <- function(de.lists, pairs, pval.field="p.value", effect.field="logFC", 
    pval.type=c("any", "some", "all"), log.p.in=FALSE, log.p.out=log.p.in, 
    output.field=NULL, full.stats=FALSE, sorted=TRUE)
{
    if (length(de.lists)!=nrow(pairs)) {
        stop("'nrow(pairs)' must be equal to 'length(de.lists)'")
    }
    if (is.null(output.field)) {
        output.field <- if (full.stats) "stats" else effect.field
    }

    pval.type <- match.arg(pval.type)
    method <- switch(pval.type, any="simes", some="holm-middle", all="berger")

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
        } else if (!is.null(effect.field)) {
            all.effects <- lapply(cur.stats, "[[", i=effect.field)
            stat.df <- DataFrame(all.effects)
        } else {
            stat.df <- preamble[,0] # output.field is NULL, so colnames is automatically empty.
        }
        colnames(stat.df) <- sprintf("%s.%s", output.field, targets)

        # Producing the output object.
        marker.set <- cbind(preamble, stat.df)
        if (sorted) {
            marker.set <- marker.set[gene.order,,drop=FALSE]
        }
        output[[host]] <- marker.set
    }

    SimpleList(output)
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

