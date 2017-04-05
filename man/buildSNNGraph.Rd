\name{buildSNNGraph}
\alias{buildSNNGraph}
\alias{buildSNNGraph,matrix-method}
\alias{buildSNNGraph,SCESet-method}

\title{Build a SNN graph}
\description{Build a shared-nearest-neighbors graph for cells based on their expression profiles.}

\usage{
\S4method{buildSNNGraph}{matrix}(x, k=10, subset.row=NULL)

\S4method{buildSNNGraph}{SCESet}(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{x}{A SCESet object, or a matrix containing expression values for each gene (row) in each cell (column).}
\item{k}{An integer scalar specifying the number of neighbors to consider during graph construction.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
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
This is faster/more memory efficient than hierarchical clustering for large numbers of cells, as it avoids the need to construct a distance matrix.
The choice of \code{k} can be roughly interpreted as the minimum cluster size.
Note that the generation of the graph uses \code{\link[FNN]{get.knn}} and may take some time.

Expression values should typically be on the log-scale, e.g., log-transformed counts.
Ranks can also be used for greater robustness, e.g., from \code{\link{quickCluster}} with \code{get.ranks=TRUE}.

Any setting of \code{subset.row} will override \code{get.spikes} in \code{buildSNNGraph,SCESet-method}.

Note that the setting of \code{k} here is slightly different from that used in SNN-Cliq.
The original implementation considers each cell to be its first nearest neighbor that contributes to \code{k}.
In \code{buildSNNGraph}, the \code{k} nearest neighbours refers to the number of \emph{other} cells.
}

\value{
An igraph graph, where nodes are cells and weighted edges represent connections between nearest neighbors.
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
