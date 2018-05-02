\name{Quick clustering}
\alias{quickCluster}
\alias{quickCluster,ANY-method}
\alias{quickCluster,SingleCellExperiment-method}

\title{Quick clustering of cells}
\description{Cluster similar cells based on rank correlations in their gene expression profiles.}

\usage{
\S4method{quickCluster}{ANY}(x, min.size=200, max.size=NULL, method=c("hclust", "igraph"),
    pc.approx=TRUE, get.ranks=FALSE, subset.row=NULL, min.mean=1, ...)

\S4method{quickCluster}{SingleCellExperiment}(x, subset.row=NULL, ..., assay.type="counts", get.spikes=FALSE)
}

\arguments{
\item{x}{
    A numeric count matrix where rows are genes and columns are cells.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{min.size}{An integer scalar specifying the minimum size of each cluster.}
\item{max.size}{An integer scalar specifying the maximum size of each cluster.}
\item{get.ranks}{A logical scalar specifying whether a matrix of adjusted ranks should be returned.}
\item{method}{A string specifying the clustering method to use.}
\item{pc.approx}{Argument passed to \code{\link{buildSNNGraph}} when \code{method="igraph"}, otherwise ignored.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{min.mean}{A numeric scalar specifying the filter to be applied on the average count for each filter prior to computing ranks.
Disabled by setting to \code{NULL}.}
\item{...}{
    For \code{quickCluster,ANY-method}, additional arguments to be passed to \code{\link{cutreeDynamic}} for \code{method="hclust"}, 
        or \code{\link{buildSNNGraph}} for \code{method="igraph"}.
    For \code{quickCluster,SingleCellExperiment-method}, additional arguments to pass to \code{quickCluster,ANY-method}.
}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
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
This should be used in situations where there are too many cells to build a distance matrix.
}
A correlation-based approach is preferred here as it is invariant to scaling normalization.
This avoids circularity between normalization and clustering, e.g., in \code{\link{computeSumFactors}}.
}

\section{Enforcing cluster sizes}{
With \code{method="hclust"}, \code{\link{cutreeDynamic}} is used to ensure that all clusters contain a minimum number of cells.
However, some cells may not be assigned to any cluster and are assigned identities of \code{"0"} in the output vector.
In most cases, this is because those cells belong in a separate cluster with fewer than \code{min.size} cells.
The function will not be able to call this as a cluster as the minimum threshold on the number of cells has not been passed.
Users are advised to check that the unassigned cells do indeed form their own cluster.
Otherwise, it may be necessary to use a different clustering algorithm.

When using \code{method="igraph"}, clusters are first identified using \code{\link[igraph]{cluster_fast_greedy}}.
If the smallest cluster contains fewer cells than \code{min.size}, it is merged with the closest neighbouring cluster.
In particular, the function will attempt to merge the smallest cluster with each other cluster.
The merge that maximizes the modularity score is selected, and a new merged cluster is formed.
This process is repeated until all (merged) clusters are larger than \code{min.size}.
}

\section{Gene selection}{
In \code{quickCluster,SingleCellExperiment-method}, spike-in transcripts are not used by default as they provide little information on the biological similarities between cells.
This may not be the case if subpopulations differ by total RNA content, in which case setting \code{get.spikes=TRUE} may provide more discriminative power.

Users can also set \code{subset.row} to specify which rows of \code{x} are to be used to calculate correlations.
This is equivalent to but more efficient than subsetting \code{x} directly, as it avoids constructing a (potentially large) temporary matrix.
Note that if \code{subset.row} is specified, it will intersect with any setting of \code{get.spikes}.

By default, the function will also filter out genes with average counts (as defined by \code{\link{calcAverage}}) below \code{min.mean}.
This removes low-abundance genes with many tied ranks, especially due to zeros, which may reduce the precision of the clustering.
We suggest setting \code{min.mean} to 1 for read count data and 0.1 for UMI data.
This can be disabled completely by setting it to \code{NULL}.
}

\section{Obtaining the scaled and centred ranks}{
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
