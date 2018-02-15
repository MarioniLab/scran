\name{findMarkers}
\alias{findMarkers}
\alias{findMarkers,ANY-method}
\alias{findMarkers,SingleCellExperiment-method}

\title{Find marker genes}
\description{Find candidate marker genes for clusters of cells, by testing for differential expression between clusters.}

\usage{
\S4method{findMarkers}{ANY}(x, clusters, block=NULL, design=NULL,
    pval.type=c("any", "all"), direction=c("any", "up", "down"), 
    lfc=0, log.p=FALSE, full.stats=FALSE, subset.row=NULL)

\S4method{findMarkers}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts", 
    get.spikes=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix-like object of normalized log-expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{clusters}{
A vector of cluster identities for all cells.
}
\item{block}{A factor specifying the blocking level for each cell.}
\item{design}{A numeric matrix containing blocking terms, i.e., uninteresting factors driving expression across cells.}
\item{pval.type}{A string specifying the type of combined p-value to be computed, i.e., Simes' or IUT.}
\item{direction}{A string specifying the direction of log-fold changes to be considered for each cluster.}
\item{lfc}{A positive numeric scalar specifying the log-fold change threshold to be tested against.}
\item{log.p}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{full.stats}{A logical scalar indicating whether all statistics (i.e., raw and BH-adjusted p-values) should be returned for each pairwise comparison.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to the matrix method.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\details{
This function performs t-tests to identify differentially expressed genes (DEGs) between pairs of clusters.
For each cluster, the log-fold changes and other statistics from all relevant pairwise comparisons are combined into a single table.
A list of such tables is returned for all clusters to define a set of potential marker genes.

Users can specify the genes to check for differential expression (DE) by setting the \code{subset.row} argument.
In addition, spike-in transcripts are ignored in the SingleCellExperiment method when \code{get.spikes=FALSE}.
If this is set, it will intersect with any non-\code{NULL} value of \code{subset.row}, i.e., only non-spike-in transcripts in the specified set will be used.
}

\section{Explanation of the hypothesis tests}{
By default, this function will perform a Welch t-test to identify DEGs between each pair of clusters.
If one of the clusters contains only one cell, a Student's t-test will be performed instead.
If both of the clusters contain only one cell, no p-value will be reported for use in ranking (see below).

If \code{block} is specified, the same t-tests are performed between clusters within each level of \code{block}.
For each pair of clusters, the p-values for each gene across all levels of \code{block} are combined using Stouffer's Z-score method.
The p-value for each level is assigned a weight inversely proportional to the expected variance of the log-fold change estimate for that level.
Blocking levels are ignored if no p-value was reported.

If \code{design} is specified, a linear model is instead fitted to the expression profile for each gene.
This linear model will include the \code{clusters} as well as any blocking factors in \code{design}.
A t-test is then performed to identify DEGs between pairs of clusters, using the values of the relevant coefficients and the gene-wise residual variance.

Note that \code{block} will override any \code{design} if both are specified.
This reflects our preference for the former, which is more robust to misspecification of the clusters, 
differences in the variance of expression in each cluster, and differences in the log-fold changes between clusters across blocks.
Nonetheless, use of \code{design} is unavoidable when blocking on real-valued covariates.
}

\section{Direction and magnitude of the log-fold change}{ 
If \code{direction="any"}, two-sided tests will be performed for each pairwise comparisons between clusters.
Otherwise, one-sided tests in the specified direction will be used to compute p-values for each gene.
This can be used to focus on genes that are upregulated in each cluster of interest, which is often easier to interpret.

To interpret the setting of \code{direction}, consider the DataFrame for cluster X, in which we are comparing to another cluster Y.
If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in cluster X compared to Y.
If \code{direction="down"}, genes will only be significant if they are downregulated in cluster X compared to Y.

The magnitude of the log-fold changes can also be tested by setting \code{lfc}.
By default, \code{lfc=0} meaning that we will reject the null upon detecting any differential expression.
If this is set to some other positive value, the null hypothesis will change depending on \code{direction}:
\itemize{
\item If \code{direction="any"}, the null hypothesis is that the true log-fold change is either \code{-lfc} or \code{lfc} with equal probability.
A two-sided p-value is computed against this composite null.
\item If \code{direction="up"}, the null hypothesis is that the true log-fold change is \code{lfc}, and a one-sided p-value is computed.
\item If \code{direction="down"}, the null hypothesis is that the true log-fold change is \code{-lfc}, and a one-sided p-value is computed.
}
This is similar to the approach used in \code{\link[limma]{treat}} and allows users to focus on genes with strong log-fold changes.
}

\section{Consolidating p-values into a ranking}{
By default, each table is sorted by the \code{Top} value when \code{pval.type="any"}.
This is the minimum rank across all pairwise comparisons for each gene, and specifies the size of the candidate marker set.
Taking all rows with \code{Top} values no greater than some integer X will yield a set containing the top X genes (ranked by significance) from each pairwise comparison.
For example, if X is 5, the set will consist of the \emph{union} of the top 5 genes from each pairwise comparison.
The marker set for each cluster allows it to be distinguished from each other cluster based on the expression of at least one gene.

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

\section{Weighting across blocking levels}{
When \code{block} is specified, the weight for the p-value in a particular level is defined as \eqn{(1/Nx + 1/Ny)^{-1}}, 
where \eqn{Nx} and \eqn{Ny} are the number of cells in clusters X and Y, respectively, for that level. 
This is inversely proportional to the expected variance of the log-fold change, provided that all clusters and blocking levels have the same variance.

In theory, a better weighting scheme would be to use the estimated standard error of the log-fold change to compute the weight.
This would be more responsive to differences in variance between blocking levels, focusing on levels with low variance and high power.
However, this is not safe in practice as genes with many zeroes can have very low standard errors, dominating the results inappropriately.

Like the p-values, the reported log-fold change for each gene is a weighted average of log-fold changes from all levels of the blocking factor. 
The weight for each log-fold change is inversely proportional to the expected variance of the log-fold change in that level.
Unlike p-values, though, this calculation will use blocking levels where both clusters contain only one cell.
}

\value{
A named list of DataFrames.
Each DataFrame corresponds to a cluster in \code{clusters}, where rows correspond to genes and are ranked by importance.
Within the DataFrame for each cluster (e.g., cluster X), there are the following fields:
\describe{
\item{\code{Top}:}{Integer, the minimum rank across all pairwise comparisons if \code{rank.type="any"}.}
\item{\code{IUT.p}:}{Numeric, the IUT p-value across all comparisons if \code{rank.type="all"}.
This is log-transformed and reported as \code{log.IUT.p} if \code{log.p=TRUE}.}
\item{\code{FDR}:}{Numeric, the BH-adjusted p-value for each gene.
This is log-transformed and reported as \code{log.FDR} if \code{log.p=TRUE}.}
\item{\code{logFC.Y}:}{Numeric for every other cluster Y in \code{clusters}, containing the log-fold change of X over Y when \code{full.stats=FALSE}.}
\item{\code{stats.Y}:}{DataFrame for every other cluster Y in \code{clusters}, returned when \code{full.stats=TRUE}.
This contains the numeric fields \code{logFC}, the log-fold change of X over Y; 
\code{p.value} or \code{log.p.value}, the (log-transformed) p-value for the pairwise comparison between X and Y;
and \code{FDR} or \code{log.FDR}, the (log-transformed) BH-adjusted p-value.
Setting \code{log.p=TRUE} will yield the log-transformed output.}
}
Genes are ranked by the \code{Top} or \code{IUT.p} column, depending on \code{rank.type}.
} 

\author{
Aaron Lun
}

\seealso{
\code{\link{normalize}}
}

\references{
Simes RJ (1986). 
An improved Bonferroni procedure for multiple tests of significance. 
\emph{Biometrika} 73:751-754.

Berger RL and Hsu JC (1996). 
Bioequivalence trials, intersection-union tests and equivalence confidence sets.
\emph{Statist. Sci.} 11, 283-319.

Whitlock MC (2005). 
Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
\emph{J. Evol. Biol.} 18, 5:1368-73.
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
kout <- kmeans(t(logcounts(y2)), centers=2) # Any clustering method is okay.
out <- findMarkers(y2, clusters=kout$cluster)
}
