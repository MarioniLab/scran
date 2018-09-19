\name{combineMarkers}
\alias{combineMarkers}

\title{Combine DE results to a marker set}
\description{Combine pairwise comparisons between groups or clusters into a single marker set for each cluster.}

\usage{
combineMarkers(de.lists, pairs, pval.field="p.value", effect.field="logFC", 
    pval.type=c("any", "all"), log.p.in=FALSE, log.p.out=log.p.in, 
    output.field=NULL, full.stats=FALSE)
}

\arguments{
\item{de.lists}{
A list-like object where each element is a data.frame or \linkS4class{DataFrame}.
Each element should represent the results of a pairwise comparison between two groups/clusters,
in which each row should contain the statistics for a single gene/feature.
Rows should be named by the feature name in the same order for all elements.
}
\item{pairs}{
A matrix, data.frame or \linkS4class{DataFrame} with two columns and number of rows equal to the length of \code{de.lists}.
Each row should specify the pair of clusters being compared for the corresponding element of \code{de.lists}.
}
\item{pval.field}{A string specifying the column name of each element of \code{de.lists} that contains the p-value.}
\item{effect.field}{A string specifying the column name of each element of \code{de.lists} that contains the effect size.}
\item{pval.type}{A string specifying the type of combined p-value to be computed, i.e., Simes' (\code{"any"}) or IUT (\code{"all"}).}
\item{log.p.in}{A logical scalar indicating if the p-values in \code{de.lists} were log-transformed.}
\item{log.p.out}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{output.field}{A string specifying the prefix of the field names containing the effect sizes.}
\item{full.stats}{A logical scalar indicating whether all statistics in \code{de.lists} should be stored in the output for each pairwise comparison.}
}

\value{
A named \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains the consolidated results for the cluster of the same name.

Within each DataFrame (say, the DataFrame for cluster X), rows correspond to genes with the fields:
\describe{
\item{\code{Top}:}{Integer, the minimum rank across all pairwise comparisons.
Only reported if \code{rank.type="any"}.}
\item{\code{p.value}:}{Numeric, the p-value across all comparisons if \code{log.p.out=FALSE}.
This is a Simes' p-value if \code{rank.type="any"}, otherwise it is an IUT p-value.}
\item{\code{log.p.value}:}{Numeric, the log-transformed version of \code{p.value} if \code{log.p.out=TRUE}.}
\item{\code{FDR}:}{Numeric, the BH-adjusted p-value for each gene if \code{log.p.out=FALSE}.}
\item{\code{log.FDR}:}{Numeric, the log-transformed adjusted p-value for each gene if \code{log.p.out=TRUE}.}
\item{\code{logFC.Y}:}{Numeric for every other cluster Y in \code{clusters}, containing the effect size of the comparison of X to Y when \code{full.stats=FALSE}.
Note that this is named according to the \code{output.field} prefix and so may not necessarily contain the log-fold change.}
\item{\code{stats.Y}:}{DataFrame for every other cluster Y in \code{clusters}, returned when \code{full.stats=TRUE}.
This contains the same fields in the corresponding entry of \code{de.lists} for the X versus Y comparison.}
}
Genes are ranked by the \code{Top} column (if available) and then the \code{p.value} column.

DataFrames are sorted according to the order of cluster IDs in \code{pairs[,1]}.
The \code{logFC.Y} columns are sorted according to the order of cluster IDs in \code{pairs[,2]} within the corresponding level of the first cluster.

Note that DataFrames are only created for clusters present in \code{pairs[,1]}.
Clusters unique to \code{pairs[,2]} will only be present within each DataFrame as Y.
}

\details{
An obvious strategy to characterizing differences between clusters is to look for genes that are differentially expressed (DE) between them.
However, this entails a number of comparisons between all pairs of clusters to comprehensively identify genes that define each cluster.
For all pairwise comparisons involving a single cluster, we would like to consolidate the DE results into a single list of candidate marker genes.
This is the intention of the \code{combineMarkers} function.
DE statistics from any testing regime can be supplied to this function - see the Examples for how this is done with t-tests from \code{\link{pairwiseTTests}}.
}

\section{Consolidating p-values into a ranking}{
By default, each table is sorted by the \code{Top} value when \code{pval.type="any"}.
This is the minimum rank across all pairwise comparisons for each gene, and specifies the size of the candidate marker set.
Taking all rows with \code{Top} values no greater than some integer X will yield a set containing the top X genes (ranked by significance) from each pairwise comparison.
For example, if X is 5, the set will consist of the \emph{union} of the top 5 genes from each pairwise comparison.
The marker set for each cluster allows it to be distinguished from every other cluster based on the expression of at least one gene.

This approach does not explicitly favour genes that are uniquely expressed in a cluster.
Such a strategy is often too stringent, especially in cases involving overclustering or cell types defined by combinatorial gene expression.
However, if \code{pval.type="all"}, the null hypothesis is that the gene is not DE in all contrasts, and the IUT p-value is computed for each gene.
This yields a \code{IUT.p} field instead of a \code{Top} field in the output table.
Ranking based on the IUT p-value will focus on genes that are uniquely DE in that cluster.
}

\section{Correcting for multiple testing}{
When \code{pval.type="any"}, a combined p-value is calculated by consolidating p-values across contrasts for each gene using Simes' method.
This represents the evidence against the null hypothesis is that the gene is not DE in any of the contrasts.
The BH method is then applied on the combined p-values across all genes to obtain the \code{FDR} field.
The same procedure is done with \code{pval.type="all"}, but using the IUT p-values across genes instead.

If \code{log.p=TRUE}, log-transformed p-values and FDRs will be reported.
This may be useful in over-powered studies with many cells, where directly reporting the raw p-values would result in many zeroes due to the limits of machine precision.

Note that the reported FDRs are intended only as a rough measure of significance.
Properly correcting for multiple testing is not generally possible when \code{clusters} is determined from the same \code{x} used for DE testing.
}

\seealso{
\code{\link{pairwiseTTests}} and \code{\link{pairwiseWilcox}} for functions that can generate \code{de.lists} and \code{pairs}.

\code{\link{findMarkers}} and \code{\link{overlapExprs}} provide wrappers that use \code{combineMarkers} on the t-test or Wilcoxon test results.
}

\references{
Simes RJ (1986). 
An improved Bonferroni procedure for multiple tests of significance. 
\emph{Biometrika} 73:751-754.

Berger RL and Hsu JC (1996). 
Bioequivalence trials, intersection-union tests and equivalence confidence sets.
\emph{Statist. Sci.} 11, 283-319.
}

\author{
Aaron Lun
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
kout <- kmeans(t(logcounts(y2)), centers=2) # Any clustering method is okay.

out <- pairwiseTTests(logcounts(y2), clusters=paste0("Cluster", kout$cluster))
comb <- combineMarkers(out$statistics, out$pairs)
comb[["Cluster1"]]
}
