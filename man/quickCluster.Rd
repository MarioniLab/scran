\name{Quick clustering}
\alias{quickCluster}
\alias{quickCluster,matrix-method}
\alias{quickCluster,SCESet-method}

\title{Quick clustering of cells}
\description{Cluster similar cells based on rank correlations in their gene expression profiles.}

\usage{
\S4method{quickCluster}{matrix}(x, min.size=200, subset.row=NULL, get.ranks=FALSE, 
    method=c("hclust", "igraph"), ...)

\S4method{quickCluster}{SCESet}(x, subset.row=NULL, ..., assay="counts", get.spikes=FALSE)
}

\arguments{
\item{x}{
    A numeric count matrix where rows are genes and columns are cells.
    Alternatively, a SCESet object containing such a matrix.
}
\item{min.size}{An integer scalar specifying the minimum size of each cluster for \code{method="hclust"}.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{get.ranks}{A logical scalar specifying whether a matrix of adjusted ranks should be returned.}
\item{method}{A string specifying the clustering method to use.}
\item{...}{
    For \code{quickCluster,matrix-method}, additional arguments to be passed to \code{\link{cutreeDynamic}} for \code{method="hclust"}, 
        or \code{\link{buildSNNGraph}} for \code{method="igraph"}.
    For \code{quickCluster,SCESet-method}, additional arguments to pass to \code{quickCluster,matrix-method}.
}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical specifying whether spike-in transcripts should be used.}
}

\details{
This function provides a correlation-based approach to quickly define clusters of a minimum size \code{min.size}.
Two clustering strategies are available:
\itemize{
\item If \code{method="hclust"}, a distance matrix is constructed using Spearman's rho on the counts between cells.
(Some manipulation is performed to convert Spearman's rho into a proper distance metric.)
Hierarchical clustering is performed and a dynamic tree cut is used to define clusters of cells.
\item If \code{method="igraph"}, a shared nearest neighbor graph is constructed using the \code{\link{buildSNNGraph}} function.
This is used to define clusters based on highly connected communities in the graph, using the \code{\link[igraph]{cluster_fast_greedy}} function.
Again, neighbors are identified using distances based on Spearman's rho.
}
A correlation-based approach is preferred here as it is invariant to scaling normalization.
This avoids circularity between normalization and clustering, e.g., in \code{\link{computeSumFactors}}.

When using \code{\link{cutreeDynamic}} some cells may not be assigned to any cluster, and are assigned identities of \code{"0"} in the output vector.
In most cases, this is because those cells belong in a separate cluster with fewer than \code{min.size} cells.
The function will not be able to call this as a cluster as the minimum threshold on the number of cells has not been passed.
Users are advised to check that the unassigned cells do indeed form their own cluster.
Otherwise, it may be necessary to use a custom clustering algorithm.

Using \code{method="igraph"} should be used in situations where there are too many cells for construction of a distance matrix.
It can also be used in cases with few cells, though it is less customizable than \code{method="hclust"} - for example, there are no options to control the minimum cluster size.
(When there are many cells, the \dQuote{quick} in the function's name refers to the number of commands you have to call rather than the computation speed.)

In \code{quickCluster,SCESet-method}, spike-in transcripts are not used by default as they provide little information on the biological similarities between cells.
This may not be the case if subpopulations differ by total RNA content, in which case setting \code{get.spikes=TRUE} may provide more discriminative power.
Users can also set \code{subset.row} to specify which rows of \code{x} are to be used to calculate correlations.
This is equivalent to but more efficient than subsetting \code{x} directly, as it avoids constructing a (potentially large) temporary matrix.
Note that if \code{subset.row} is specified, it will overwrite any setting of \code{get.spikes}.

Users can also set \code{get.ranks=TRUE}, in which case a matrix of ranks will be returned.
Each column contains the ranks for the expression values within a single cell after standardization and mean-centring.
Computing Euclidean distances between the rank vectors for pairs of cells will yield the same correlation-based distance as that used above.
This allows users to apply their own clustering algorithms on the ranks, which protects against outliers and is invariant to scaling (at the cost of sensitivity).
}

\value{
If \code{get.ranks=FALSE}, a character vector of cluster identities for each cell in \code{counts} is returned.

If \code{get.ranks=TRUE}, a numeric matrix is returned where each column contains ranks for the expression values in each cell.
}

\author{
Aaron Lun and Karsten Bach
}

\seealso{
\code{\link{cutreeDynamic}},
\code{\link{computeSumFactors}},
\code{\link{buildSNNGraph}}
}

\examples{
set.seed(100)
popsize <- 200
ngenes <- 1000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs, size=1), ncol=popsize, byrow=TRUE)

clusters <- quickCluster(counts, min.size=20)
clusters <- quickCluster(counts, method="igraph")
}

\keyword{normalization}

\references{
van Dongen S and Enright AJ (2012).
Metric distances derived from cosine similarity and Pearson and Spearman correlations.
\emph{arXiv} 1208.3145

Lun ATL, Bach K and Marioni JC (2016).
Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
\emph{Genome Biol.} 17:75
}
