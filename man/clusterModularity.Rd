\name{clusterModularity}
\alias{clusterModularity}

\title{Compute the cluster-wise modularity}
\description{Calculate the modularity of each cluster from a graph, based on a null model of random connections between nodes.} 

\usage{
clusterModularity(graph, clusters, get.values=FALSE) 
}

\arguments{
\item{graph}{A graph object from \pkg{igraph}, like that produced by \code{\link{buildSNNGraph}}.}
\item{clusters}{A factor specifying the cluster identity for each node.}
\item{get.values}{A logical scalar indicating whether the observed and expected edge weights should be returned.}
}

\value{
If \code{get.values=FALSE}, a symmetric numeric matrix of order equal to the number of clusters is returned.
Each entry corresponds to a pair of clusters and is proportional to the difference between the observed and expected edge weights between those clusters.

If \code{get.values=TRUE}, a list is returned containing two symmetric numeric matrices. 
The \code{observed} matrix contains the observed sum of edge weights between and within clusters,
while the \code{expected} matrix contains the expected sum of edge weights under the random model.
}

\details{
This function computes a modularity score in the same manner as that from \code{\link[igraph]{modularity}}.
The modularity is defined as the difference between the observed and expected number of edges between nodes in the same cluster.
The expected number of edges is defined by a null model where edges are randomly distributed among nodes.
The same logic applies for weighted graphs, replacing the number of edges with the summed weight of edges.

Whereas \code{\link[igraph]{modularity}} returns a modularity score for the entire graph, \code{clusterModularity} provides scores for the individual clusters.
This allows users to determine which clusters are enriched for intra-cluster edges based on their high modularity scores.
For comparison, \code{clusterModularity} also reports the modularity scores between pairs of clusters.
The sum of the diagonal elements of the output matrix should be equal to the output of \code{\link[igraph]{modularity}} (after supplying weights to the latter, if necessary).
}

\author{
Aaron Lun
}

\seealso{
\code{\link{buildSNNGraph}},
\code{\link[igraph]{modularity}}
}

\examples{
example(buildSNNGraph) # using the mocked-up graph in this example.

# Examining the modularity values directly.
out <- clusterModularity(g, clusters)
image(out)

# Alternatively, compare the ratio of observed:expected.
out <- clusterModularity(g, clusters, get.values=TRUE)
log.ratio <- log2(out$observed/out$expected + 1)
image(log.ratio)
}

