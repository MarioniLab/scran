\name{build*NNGraph}

\alias{buildSNNGraph}
\alias{buildSNNGraph,ANY-method}
\alias{buildSNNGraph,SingleCellExperiment-method}

\alias{buildKNNGraph}
\alias{buildKNNGraph,ANY-method}
\alias{buildKNNGraph,SingleCellExperiment-method}

\title{Build a nearest-neighbor graph}
\description{Build a shared or k-nearest-neighbors graph for cells based on their expression profiles.}

\usage{
\S4method{buildSNNGraph}{ANY}(x, k=10, d=50, transposed=FALSE, 
    pc.approx=FALSE, rand.seed=1000, irlba.args=list(), 
    subset.row=NULL, BPPARAM=SerialParam())

\S4method{buildSNNGraph}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts", 
    get.spikes=FALSE, use.dimred=NULL)

\S4method{buildKNNGraph}{ANY}(x, k=10, d=50, directed=FALSE, transposed=FALSE, 
    pc.approx=FALSE, rand.seed=1000, irlba.args=list(), 
    subset.row=NULL, BPPARAM=SerialParam())

\S4method{buildKNNGraph}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts", 
    get.spikes=FALSE, use.dimred=NULL)
}

\arguments{
\item{x}{A SingleCellExperiment object, or a matrix containing expression values for each gene (row) in each cell (column).
If it is matrix, it can also be transposed.}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider during graph construction.}
\item{d}{An integer scalar specifying the number of dimensions to use for the k-NN search.}
\item{directed}{A logical scalar indicating whether the output of \code{buildKNNGraph} should be a directed graph.}
\item{transposed}{A logical scalar indicating whether \code{x} is transposed (i.e., rows are cells).}
\item{pc.approx}{A logical scalar indicating whether approximate PCA should be performed.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{rand.seed}{A numeric scalar specifying the seed for approximate PCA when \code{pc.approx=TRUE}.
This can be set to \code{NA} to use the existing session seed.}
\item{irlba.args}{A named list of additional arguments to pass to \code{\link[irlba]{prcomp_irlba}} when \code{pc.approx=TRUE}.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{...}{Additional arguments to pass to \code{buildSNNGraph,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{get.spikes}{A logical scalar specifying whether spike-in transcripts should be used.}
\item{use.dimred}{A string specifying whether existing values in \code{reducedDims(x)} should be used.}
}

\details{
The \code{buildSNNGraph} method builds a shared nearest-neighbour graph using cells as nodes.
For each cell, its \code{k} nearest neighbours are identified based on Euclidean distances in their expression profiles.
An edge is drawn between all pairs of cells that share at least one neighbour.
The weight of the edge between two cells is determined by the ranking of the shared nearest neighbors.
More shared neighbors, or shared neighbors that are close to both cells, will yield larger weights.

The aim is to use the SNN graph to perform community-based clustering, using various methods in the \pkg{igraph} package.
This is faster/more memory efficient than hierarchical clustering for large numbers of cells.
In particular, it avoids the need to construct a distance matrix for all pairs of cells.
The choice of \code{k} can be roughly interpreted as the minimum cluster size.

Note that the setting of \code{k} here is slightly different from that used in SNN-Cliq.
The original implementation considers each cell to be its first nearest neighbor that contributes to \code{k}.
In \code{buildSNNGraph}, the \code{k} nearest neighbours refers to the number of \emph{other} cells.

The \code{buildKNNGraph} method builds a simpler k-nearest neighbour graph.
Cells are again nodes, and edges are drawn between each cell and its k-nearest neighbours.
No weighting of the edges is performed.
In theory, these graphs are directed as nearest neighour relationships may not be reciprocal.
However, by default, \code{directed=FALSE} such that an undirected graph is returned.
}

\section{Choice of input data}{
In practice, PCA is performed on \code{x} to obtain the first \code{d} principal components.
This is necessary in order to perform the k-NN search (done using the \code{\link[FNN]{get.knn}} function) in reasonable time.
By default, the first 50 components are chosen, which should retain most of the substructure in the data set.
If \code{d} is \code{NA} or less than the number of cells, no dimensionality reduction is performed.
If \code{pc.approx=TRUE}, \code{\link[irlba]{prcomp_irlba}} will be used to quickly obtain the first \code{d} PCs.

Expression values in \code{x} should typically be on the log-scale, e.g., log-transformed counts.
Ranks can also be used for greater robustness, e.g., from \code{\link{quickCluster}} with \code{get.ranks=TRUE}.
(Dimensionality reduction is still okay when ranks are provided - running PCA on ranks is equivalent to running MDS on the distance matrix derived from Spearman's rho.)
If the input matrix is already transposed, \code{transposed=TRUE} avoids an unnecessary internal transposition. 

By default, spike-in transcripts are removed from the expression matrix in \code{buildSNNGraph,SCESet-method}.
However, any non-\code{NULL} setting of \code{subset.row} will override \code{get.spikes}.
If \code{use.dimred} is not \code{NULL}, existing PCs are used from the specified entry of \code{reducedDims(x)}, 
and any setting of \code{d}, \code{subset.row} and \code{get.spikes} are ignored.
}

\value{
An igraph-type graph, where nodes are cells and edges represent connections between nearest neighbors.
For \code{buildSNNGraph}, these edges are weighted by the number of shared nearest neighbors.
For \code{buildKNNGraph}, edges are not weighted but may be directed if \code{directed=TRUE}.
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

clusters <- igraph::cluster_fast_greedy(g)$membership
table(clusters)
}
