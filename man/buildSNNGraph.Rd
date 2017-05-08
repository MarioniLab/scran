\name{buildSNNGraph}
\alias{buildSNNGraph}
\alias{buildSNNGraph,matrix-method}
\alias{buildSNNGraph,SCESet-method}

\title{Build a SNN graph}
\description{Build a shared-nearest-neighbors graph for cells based on their expression profiles.}

\usage{
\S4method{buildSNNGraph}{matrix}(x, k=10, d=50, subset.row=NULL, BPPARAM=SerialParam())

\S4method{buildSNNGraph}{SCESet}(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{x}{A SCESet object, or a matrix containing expression values for each gene (row) in each cell (column).}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider during graph construction.}
\item{d}{An integer scalar specifying the number of dimensions to use for the k-NN search.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{...}{Additional arguments to pass to \code{buildSNNGraph,matrix-method}.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical specifying whether spike-in transcripts should be used.}
}

\details{
This function builds a SNN graph using cells as nodes.
Each cell is connected to its \code{k} nearest neighbors, based on Euclidean distances in their expression profiles.
The weight of the edge between two cells is determined by the ranking of their shared nearest neighbors.
More shared neighbors, or shared neighbors that are close to both cells, will yield larger weights.

The aim is to use the SNN graph to perform community-based clustering, using various methods in the \pkg{igraph} package.
This is faster/more memory efficient than hierarchical clustering for large numbers of cells.
In particular, it avoids the need to construct a distance matrix for all pairs of cells.
The choice of \code{k} can be roughly interpreted as the minimum cluster size.

In practice, PCA is performed on \code{x} to obtain the first \code{d} principal components.
This is necessary in order to perform the k-NN search (done using the \code{\link[FNN]{get.knn}} function) in reasonable time.
By default, the first 50 components are chosen, which should retain most of the substructure in the data set.
If \code{d} is \code{NA} or less than the number of cells, no dimensionality reduction is performed.

Expression values should typically be on the log-scale, e.g., log-transformed counts.
Ranks can also be used for greater robustness, e.g., from \code{\link{quickCluster}} with \code{get.ranks=TRUE}.
(Dimensionality reduction is still okay when ranks are provided - running PCA on ranks is equivalent to running MDS on the distance matrix derived from Spearman's rho.)

By default, spike-in transcripts are removed from the expression matrix in \code{buildSNNGraph,SCESet-method}.
However, any non-\code{NULL} setting of \code{subset.row} will override \code{get.spikes}.

Note that the setting of \code{k} here is slightly different from that used in SNN-Cliq.
The original implementation considers each cell to be its first nearest neighbor that contributes to \code{k}.
In \code{buildSNNGraph}, the \code{k} nearest neighbours refers to the number of \emph{other} cells.
}

\value{
An igraph-type graph, where nodes are cells and weighted edges represent connections between nearest neighbors.
}

\author{
Aaron Lun
}

\seealso{
    \code{\link[FNN]{get.knn}},
    \code{\link[igraph]{make_graph}}
}

\references{
Xu C and Su Z (2015).
Identification of cell types from single-cell transcriptomes using a novel clustering method.
\emph{Bioinformatics} 31:1974-80
}

\examples{
exprs <- matrix(rnorm(100000), ncol=100)
g <- buildSNNGraph(exprs)

library(igraph) # lots of algorithms can be used
clusters <- cluster_fast_greedy(g)$membership
}
