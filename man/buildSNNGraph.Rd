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
\S4method{buildSNNGraph}{ANY}(x, k=10, d=50, type=c("rank", "number"), transposed=FALSE,
    pc.approx=NULL, irlba.args=list(), subset.row=NULL, BNPARAM=KmknnParam(), 
    BSPARAM=ExactParam(), BPPARAM=SerialParam())

\S4method{buildSNNGraph}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts", 
    get.spikes=FALSE, use.dimred=NULL)

\S4method{buildKNNGraph}{ANY}(x, k=10, d=50, directed=FALSE, transposed=FALSE, 
    pc.approx=NULL, irlba.args=list(), subset.row=NULL, BNPARAM=KmknnParam(), 
    BSPARAM=ExactParam(), BPPARAM=SerialParam())

\S4method{buildKNNGraph}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts", 
    get.spikes=FALSE, use.dimred=NULL)
}

\arguments{
\item{x}{A SingleCellExperiment object, or a matrix containing expression values for each gene (row) in each cell (column).
If it is matrix, it can also be transposed.}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider during graph construction.}
\item{d}{An integer scalar specifying the number of dimensions to use for the k-NN search.}
\item{type}{A string specifying the type of weighting scheme to use for shared neighbors.}
\item{directed}{A logical scalar indicating whether the output of \code{buildKNNGraph} should be a directed graph.}
\item{transposed}{A logical scalar indicating whether \code{x} is transposed (i.e., rows are cells).}
\item{pc.approx, irlba.args}{Deprecated, use \code{BSPARAM} instead.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{BNPARAM}{A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.}
\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, if \code{d} is not \code{NA}.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object to use for parallel processing.}
\item{...}{Additional arguments to pass to \code{buildSNNGraph,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{get.spikes}{See \code{?"\link{scran-gene-selection}"}.}
\item{use.dimred}{A string specifying whether existing values in \code{reducedDims(x)} should be used.}
}

\details{
The \code{buildSNNGraph} method builds a shared nearest-neighbour graph using cells as nodes.
For each cell, its \code{k} nearest neighbours are identified based on Euclidean distances in their expression profiles.
An edge is drawn between all pairs of cells that share at least one neighbour, weighted by the characteristics of the shared nearest neighbors:
\itemize{
    \item If \code{type="rank"}, the weighting scheme defined by Xu and Su (2015) is used.
        The weight between two nodes is equal to \eqn{k - r/2}, where \eqn{r} is the smallest sum of ranks for any shared neighboring node.
        For example, if both nodes had the same closest neighbor, the weight would be \eqn{k - 1}.
        For the purposes of this ranking, each node has a rank of zero in its own nearest-neighbor set. 
    \item If \code{type="number"}, the weight between two nodes is simply the number of shared nearest neighbors between them.
        The weight can range from zero to \eqn{k + 1}, as the node itself is included in its own nearest-neighbor set.
        This is a simpler scheme that is also slightly faster but does not account for the ranking of neighbors within each set.
}
More shared neighbors, or shared neighbors that are close to both cells, will generally yield larger weights.

The aim is to use the SNN graph to perform clustering of cells via community detection algorithms in the \pkg{igraph} package.
This is faster and more memory efficient than hierarchical clustering for large numbers of cells.
In particular, it avoids the need to construct a distance matrix for all pairs of cells.
Only the identities of nearest neighbours are required, which can be obtained quickly with methods in the \pkg{BiocNeighbors} package.

The choice of \code{k} can be roughly interpreted as the minimum cluster size.
Smaller values of \code{k} will generally yield smaller, more resolved clusters upon running community detection algorithms.
By comparison, increasing \code{k} will increase the connectivity of the graph and make it more difficult to resolve different communities.

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
This is necessary in order to perform the k-NN search (done using the \code{\link{findKNN}} function) in reasonable time.
By default, the first 50 components are chosen, which should retain most of the substructure in the data set.
If \code{d} is \code{NA} or greater than or equal to the number of cells, no dimensionality reduction is performed.
    
The PCA is performed using methods the \code{\link{runSVD}} function from the \pkg{BiocSingular} package.
To improve speed, this can be done using approximate algorithms by modifying \code{BSPARAM}, e.g., to \code{\link{IrlbaParam}()}.
Approximate algorithms will converge towards the correct result but often involve some random initialization and thus are technically dependent on the session seed.
For full reproducibility, users are advised to call \code{\link{set.seed}} beforehand when using this option.

Expression values in \code{x} should typically be on the log-scale, e.g., log-transformed counts.
Ranks can also be used for greater robustness, e.g., from \code{\link{quickCluster}} with \code{get.ranks=TRUE}.
(Dimensionality reduction is still okay when ranks are provided - running PCA on ranks is equivalent to running MDS on the distance matrix derived from Spearman's rho.)

If the input matrix \code{x} is already transposed for the ANY method, \code{transposed=TRUE} avoids an unnecessary internal transposition.
A typical use case is when \code{x} contains some reduced dimension coordinates with cells in the rows.
In such cases, setting \code{transposed=TRUE} and \code{d=NA} will use the input coordinates directly for graph-building.

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
See \code{\link{make_graph}} for details on the graph output object.

See \code{\link{cluster_walktrap}}, \code{\link{cluster_louvain}} and related functions in \pkg{igraph} for clustering based on the produced graph.

Also see \code{\link{findKNN}} for specifics of the nearest-neighbor search.
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
