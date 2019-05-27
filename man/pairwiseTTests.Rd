\name{pairwiseTTests}
\alias{pairwiseTTests}

\title{Perform pairwise t-tests}
\description{Perform pairwise Welch t-tests between groups of cells, possibly after blocking on uninteresting factors of variation.}

\usage{
pairwiseTTests(x, clusters, block=NULL, design=NULL, 
	direction=c("any", "up", "down"), lfc=0, log.p=FALSE, 
	gene.names=rownames(x), subset.row=NULL, BPPARAM=SerialParam())
}

\arguments{
\item{x}{A numeric matrix-like object of normalized log-expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.}
\item{clusters}{A vector of cluster identities for all cells.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{design}{A numeric matrix containing blocking terms for uninteresting factors.
Note that these should not be confounded with \code{clusters} or contain an intercept, see Details.}
\item{direction}{A string specifying the direction of log-fold changes to be considered for each cluster.}
\item{lfc}{A positive numeric scalar specifying the log-fold change threshold to be tested against.}
\item{log.p}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{gene.names}{A character vector of gene names with one value for each row of \code{x}.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{BPPARAM}{A BiocParallelParam object indicating whether and how parallelization should be performed across genes.}
}

\details{
This function performs t-tests to identify differentially expressed genes (DEGs) between pairs of clusters.
A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of clusters.
This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.

By default, this function will perform a Welch t-test to identify DEGs between each pair of clusters.
This is simple, fast and performs reasonably well for single-cell count data (Soneson and Robinson, 2018).
However, if one of the clusters contains fewer than two cells, no p-value will be reported for comparisons involving that cluster.
A warning will also be raised about insufficient degrees of freedom (d.f.) in such cases.

If \code{block} is specified, t-tests are performed between clusters within each level of \code{block}.
For each pair of clusters, the p-values for each gene across all levels of \code{block} are combined using Stouffer's Z-score method.
The p-value for each level is assigned a weight inversely proportional to the expected variance of the log-fold change estimate for that level.
Blocking levels are ignored if no p-value was reported, e.g., if there were insufficient cells for a cluster in a particular level. 
Comparisons may also yield \code{NA} p-values (along with a warning about the lack of d.f.) if the two clusters do not co-occur in the same block.

If \code{design} is specified, a linear model is instead fitted to the expression profile for each gene.
This linear model will include the \code{clusters} as well as any blocking factors in \code{design}.
A t-test is then performed to identify DEGs between pairs of clusters, using the values of the relevant coefficients and the gene-wise residual variance.
Note that \code{design} must be full rank when combined with the \code{clusters} terms, i.e., there should not be any confounding variables.
Similarly, any intercept column should be removed beforehand.

Note that \code{block} will override any \code{design} if both are specified.
This reflects our preference for the former, which accommodates differences in the variance of expression in each cluster via Welch's t-test.
As a result, it is more robust to misspecification of the clusters, as misspecified clusters (and inflated variances) do not affect the inferences for other clusters.
Use of \code{block} also avoids assuming additivity of effects between the blocking factors and the cluster identities.

Nonetheless, use of \code{design} is unavoidable when blocking on real-valued covariates.
It is also useful for ensuring that log-fold changes/p-values are computed for comparisons between all pairs of clusters
(assuming that \code{design} is not confounded with the cluster identities).
This may not be the case with \code{block} if a pair of clusters never co-occur in a single blocking level. 
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
This is similar to the approach used in \code{\link[limma:eBayes]{treat}} and allows users to focus on genes with strong log-fold changes.
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
A list is returned containing \code{statistics} and \code{pairs}.

The \code{statistics} element is itself a list of \linkS4class{DataFrame}s.
Each DataFrame contains the statistics for a comparison between a pair of clusters,
including the log-fold changes, p-values and false discovery rates.

The \code{pairs} element is a DataFrame with one row corresponding to each entry of \code{statistics}.
This contains the fields \code{first} and \code{second}, 
specifying the two clusters under comparison in the corresponding DataFrame in \code{statistics}.

In each DataFrame in \code{statistics}, the log-fold change represents the change in the \code{first} cluster compared to the \code{second} cluster.
Note that switching the \code{first} and \code{second} clusters will affect the sign of the log-fold change and, when \code{direction!="any"}, the size of the p-value itself.
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

Lun ATL (2018).
Comments on marker detection in \emph{scran}.
\url{https://ltla.github.io/SingleCellThoughts/software/marker_detection/comments.html}
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
kout <- kmeans(t(logcounts(y2)), centers=2) # Any clustering method is okay.

# Vanilla application:
out <- pairwiseTTests(logcounts(y2), clusters=kout$cluster)
out

# Directional with log-fold change threshold:
out <- pairwiseTTests(logcounts(y2), clusters=kout$cluster, direction="up", lfc=0.2)
out
}

