\name{doubletCluster}
\alias{doubletCluster}

\title{Detect doublet clusters}
\description{Identify potential clusters of doublet cells based on intermediate expression profiles.}

\usage{
\S4method{doubletCluster}{ANY}(x, clusters, subset.row=NULL, threshold=0.05, ...) 

\S4method{doubletCluster}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="counts", 
    get.spikes=FALSE)
}

\arguments{
\item{x}{
A numeric matrix-like object of count values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{clusters}{
    A vector of cluster identities for all cells.
}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{threshold}{A numeric scalar specifying the FDR threshold with which to identify significant genes.}
\item{...}{For the ANY method, additional arguments to pass to \code{\link{findMarkers}}. 

For the SingleCellExperiment method, additional arguments to pass to the ANY method.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\value{
A DataFrame containing one row per cluster with the following fields:
\describe{
    \item{\code{N}:}{Integer, number of genes that are significantly non-intermediate in the current cluster compared to the two putative source clusters.}
    \item{\code{source1}:}{String specifying the identity of the first source cluster.}
    \item{\code{source2}:}{String specifying the identity of the second source cluster.}
    \item{\code{best}:}{String specifying the identify of the top gene with the lowest p-value against the doublet hypothesis for this combination of current and source clusters.}
    \item{\code{p.value}:}{Numeric, containing the adjusted p-value for the \code{best} gene.} 
    \item{\code{prop}:}{Numeric, proportion of cells in the current cluster.}
    \item{\code{lib.size1}:}{Numeric, ratio of the median library sizes for the first source cluster to the current cluster.}
    \item{\code{lib.size2}:}{Numeric, ratio of the median library sizes for the second source cluster to the current cluster.}
}
Each row is named according to its cluster.
}

\details{
For each cluster, we examine all possible pairs of \dQuote{source} clusters, hypothesizing that the current cluster is a doublet formed from the two sources.
If so, gene expression in the current cluster should be strictly intermediate between the two sources after library size normalization.
This approach is similar to the method used by Bach et al. to detect clusters of doublet cells.

We apply pairwise t-tests to the normalized log-expression profiles (see \code{\link[scater]{normalize}}) to reject this null hypothesis.
This is done by identifying genes that are consistently up- or down-regulated in the current cluster compared to \emph{both} of the sources.
We count the number of genes that reject the null hypothesis at the specified FDR \code{threshold}.
For each cluster, the most likely pair of source clusters is that which minimizes the number of significant genes.

Potential doublet clusters are identified using the following characteristics:
\itemize{
    \item Low number of significant genes, i.e., \code{N} in the output DataFrame.
The threshold can be identified by looking for small outliers in \code{log(N)} across all clusters,
under the assumption that most clusters are \emph{not} doublets (and thus should have high \code{N}).
    \item A reasonable proportion of cells in the cluster, i.e., \code{prop}.
This requires some expectation of the doublet rate in the experimental protocol. 
    \item Library sizes of the source clusters that are below that of the current cluster, i.e., \code{lib.size*} values below unity.
This assumes that the doublet cluster will contain more RNA and have more counts than either of the two source clusters.    
}

The reported \code{p.value} is of little use in a statistical sense, and is only provided for inspection.
Technically, it could be treated as the Simes combined p-value against the doublet hypothesis for the current cluster.
However, this does not account for the multiple testing across all pairs of clusters for each chosen cluster, 
especially as we are chosing the pair that is most concordant with the doublet null hypothesis.
}

\author{
Aaron Lun
}

\references{
Bach K, Pensa S, Grzelak M, Hadfield J, Adams DJ, Marioni JC and Khaled WT (2017). 
Differentiation dynamics of mammary epithelial cells revealed by single-cell RNA sequencing. 
\emph{Nat Commun.} 8, 1:2128.
}

\examples{
# Mocking up an example.
ngenes <- 100
mu1 <- 2^rexp(ngenes)
mu2 <- 2^rnorm(ngenes)

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes)
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes)
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes)

counts <- cbind(counts.1, counts.2, counts.m)
clusters <- rep(1:3, c(ncol(counts.1), ncol(counts.2), ncol(counts.m)))

# Find potential doublets.
dbl <- doubletCluster(counts, clusters)

library(scater)
possible.doublet <- isOutlier(dbl$N, log=TRUE, type="lower")
}
