\name{pairwiseWilcox}
\alias{pairwiseWilcox}

\title{Perform pairwise Wilcoxon rank sum tests}
\description{Perform pairwise Wilcoxon rank sum tests between groups of cells, possibly after blocking on uninteresting factors of variation.}

\usage{
pairwiseWilcox(x, clusters, block=NULL, direction=c("any", "up", "down"),
    log.p=FALSE, gene.names=rownames(x), subset.row=NULL, tol=1e-8,
    BPPARAM=SerialParam())
}

\arguments{
\item{x}{A numeric matrix-like object of normalized log-expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.}
\item{clusters}{A vector of cluster identities for all cells.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{direction}{A string specifying the direction of effects to be considered for each cluster.}
\item{log.p}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{gene.names}{A character vector of gene names with one value for each row of \code{x}.}
\item{subset.row}{A logical, integer or character scalar specifying the rows of \code{x} to use.}
\item{tol}{Numeric scalar specifying the tolerance for tied values when \code{x} is numeric.}
\item{BPPARAM}{A BiocParallelParam object indicating whether and how parallelization should be performed across genes.}
}

\details{
This function performs Wilcoxon rank sum tests to identify differentially expressed genes (DEGs) between pairs of clusters.
A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of clusters.
This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.
Users can specify the genes to check for differential expression by setting the \code{subset.row} argument.

Effect sizes represent the probability of sampling an observation in one cluster that is \emph{greater} than a random value sampled from the other cluster.
Thus, effect sizes close to 1 indicate that the probability mass of the first cluster is completely above that of the second cluster.
Conversely, effect sizes close to 0 mean that all observations in the first cluster are below those of the second cluster.

Wilcoxon rank sum tests are more robust to outliers and insensitive to non-normality, in contrast to t-tests in \code{\link{pairwiseTTests}}.
However, they take longer to run, the effect sizes are less interpretable, and there are more subtle violations of its assumptions in real data.
For example, the i.i.d. assumptions are unlikely to hold after scaling normalization due to differences in variance.
}

\section{Blocking on uninteresting factors}{
If \code{block} is specified, the same t-tests are performed between clusters within each level of \code{block}.
For each pair of clusters, the p-values for each gene across all levels of \code{block} are combined using Stouffer's Z-score method.
The p-value for each level is assigned a weight inversely proportional to the expected variance of the log-fold change estimate for that level.
Blocking levels are ignored if no p-value was reported, e.g., if there were insufficient cells for a cluster in a particular level. 

The weight for the p-value in a particular level of \code{block} is defined as \eqn{N_xN_y},
where \eqn{N_x} and \eqn{N_y} are the number of cells in clusters X and Y, respectively, for that level. 
This means that p-values from blocks with more cells will have a greater contribution to the combined p-value for each gene.
}

\section{Direction of the effect}{
If \code{direction="any"}, two-sided Wilcoxon rank sum tests will be performed for each pairwise comparisons between clusters.
Otherwise, one-sided tests in the specified direction will be used to compute p-values for each gene.
This can be used to focus on genes that are upregulated in each cluster of interest, which is often easier to interpret.

To interpret the setting of \code{direction}, consider the DataFrame for cluster X, in which we are comparing to another cluster Y.
If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in cluster X compared to Y.
If \code{direction="down"}, genes will only be significant if they are downregulated in cluster X compared to Y.
}

\value{
A list is returned containing \code{statistics}, a list of \linkS4class{DataFrame}s where each element contains the statistics for a comparison between a pair of clusters;
and \code{pairs}, a DataFrame specifying the identity of pairs being compared in each entry of \code{statistics}.

In each DataFrame, the effect size represents the probability of sampling a value in the \code{first} cluster greater than a random value from the \code{second} cluster.
Note that the order of clusters will affect the sign of the log-fold change and, when \code{direction!="any"}, the size of the p-value itself.
}

\author{
Aaron Lun
}

\references{
Whitlock MC (2005). 
Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
\emph{J. Evol. Biol.} 18, 5:1368-73.

Soneson C and Robinson MD (2018). 
Bias, robustness and scalability in single-cell differential expression analysis. 
\emph{Nat. Methods}
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
kout <- kmeans(t(logcounts(y2)), centers=2) # Any clustering method is okay.

# Vanilla application:
out <- pairwiseWilcox(logcounts(y2), clusters=kout$cluster)
out

# Directional:
out <- pairwiseWilcox(logcounts(y2), clusters=kout$cluster, direction="up")
out
}

