\name{pairwiseBinom}
\alias{pairwiseBinom}

\title{Perform pairwise binomial tests}
\description{Perform pairwise binomial tests between groups of cells, 
possibly after blocking on uninteresting factors of variation.}

\usage{
pairwiseBinom(x, clusters, block=NULL, direction=c("any", "up", "down"),
    log.p=FALSE, gene.names=rownames(x), subset.row=NULL, threshold=1e-8,
    BPPARAM=SerialParam())
}

\arguments{
\item{x}{A numeric matrix-like object of counts,
where each column corresponds to a cell and each row corresponds to a gene.}
\item{clusters}{A vector of cluster identities for all cells.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{direction}{A string specifying the direction of effects to be considered for each cluster.}
\item{log.p}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{gene.names}{A character vector of gene names with one value for each row of \code{x}.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{threshold}{Numeric scalar specifying the value below which a gene is presumed to be not expressed.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object indicating how parallelization should be performed across genes.}
}

\details{
This function performs exact binomial tests to identify marker genes between pairs of clusters.
Here, the null hypothesis is that the proportion of cells expressing a gene is the same between clusters.
A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of clusters.
This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.

Effect sizes for each comparison are reported as log2-odds ratios,
i.e., the ratio of the proportion of expressing cells in one cluster to the proportion in another cluster.
Large log-ORs correspond to large relative differences in these proportions,
where the sign indicates the direction of the change in proportions.
We add a pseudo-count of 0.5 to all counts prior to computing the log-ORs to avoid undefined values.

\code{x} can be a count matrix or any transformed counterpart where zeroes remain zero and non-zeroes remain non-zero.
This is true of any scaling normalization and monotonic transformation like the log-transform.
If the transformation breaks this rule, some adjustment of \code{threshold} is necessary.

A consequence of the transformation-agnostic behaviour of this function is that it will not respond to normalization.
Differences in library size will not be considered by this function.
However, this is not necessarily problematic for marker gene detection -
users can treat this as \emph{retaining} information about the total RNA content, analogous to spike-in normalization.
}

\section{Blocking on uninteresting factors}{
If \code{block} is specified, binomial tests are performed between clusters within each level of \code{block}.
For each pair of clusters, the p-values for each gene across 
all levels of \code{block} are combined using Stouffer's weighted Z-score method.

The weight for the p-value in a particular level of \code{block} is defined as \eqn{N_xN_y},
where \eqn{N_x} and \eqn{N_y} are the number of cells in clusters X and Y, respectively, for that level. 
This means that p-values from blocks with more cells will have a greater contribution to the combined p-value for each gene.

Blocking levels are ignored if no p-value was reported for a particular pair of clusters, 
e.g., if there were no cells for either cluster in a particular level. 
}

\section{Direction of the effect}{
If \code{direction="any"}, approximate two-sided binomial tests will be performed for each pairwise comparisons between clusters.
For other \code{direction}, one-sided tests in the specified direction will be used to compute p-values for each gene.
This can be used to focus on genes that are upregulated in each cluster of interest, which is often easier to interpret.

The approximation of a two-sided test is done by combining two one-sided tests using a Bonferroni correction.
This is done for various logistical purposes;
it is also the only way to combine p-values across blocks in a direction-aware manner.
As a result, the two-sided p-value reported here will not be the same as that from \code{\link{binom.test}}.
In practice, they are usually similar enough that this is not a major concern.

To interpret the setting of \code{direction}, consider the DataFrame for cluster X, in which we are comparing to another cluster Y.
If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in cluster X compared to Y.
If \code{direction="down"}, genes will only be significant if they are downregulated in cluster X compared to Y.
See \code{?\link{binom.test}} for more details on the interpretation of one-sided Wilcoxon rank sum tests.

When combining across batches, one-sided p-values in the same direction are combined first.
Then, if \code{direction="any"}, the two combined p-values from both directions are combined.
This ensures that a gene only receives a low overall p-value if it changes in the same direction across batches.
}

\value{
A list is returned containing \code{statistics} and \code{pairs}.

The \code{statistics} element is itself a list of \linkS4class{DataFrame}s.
Each DataFrame contains the statistics for a comparison between a pair of clusters,
including the overlap proportions, p-values and false discovery rates.

The \code{pairs} element is a DataFrame with one row corresponding to each entry of \code{statistics}.
This contains the fields \code{first} and \code{second}, 
specifying the two clusters under comparison in the corresponding DataFrame in \code{statistics}.

In each DataFrame in \code{statistics}, the log-odds ratio represents the log-ration of the proportion of expressing cells in the \code{first} cluster compared to the expressing proportion in the \code{second} cluster.
Switching the \code{first} and \code{second} clusters will affect the value of the overlap and, when \code{direction!="any"}, the size of the p-value itself.
}

\author{
Aaron Lun
}

\references{
Whitlock MC (2005). 
Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
\emph{J. Evol. Biol.} 18, 5:1368-73.
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
kout <- kmeans(t(logcounts(y2)), centers=2) # Any clustering method is okay.

# Vanilla application:
out <- pairwiseBinom(logcounts(y2), clusters=kout$cluster)
out

# Directional:
out <- pairwiseBinom(logcounts(y2), clusters=kout$cluster, direction="up")
out
}

