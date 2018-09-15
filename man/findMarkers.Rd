\name{findMarkers}
\alias{findMarkers}
\alias{findMarkers,ANY-method}
\alias{findMarkers,SingleCellExperiment-method}

\title{Find marker genes}
\description{Find candidate marker genes for clusters of cells, by testing for differential expression between clusters.}

\usage{
\S4method{findMarkers}{ANY}(x, clusters, gene.names=rownames(x), block=NULL, design=NULL,
    pval.type=c("any", "all"), direction=c("any", "up", "down"), lfc=0, 
    log.p=FALSE, full.stats=FALSE, subset.row=NULL)

\S4method{findMarkers}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts", 
    get.spikes=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix-like object of normalized log-expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{clusters}{A vector of cluster identities for all cells.}
\item{gene.names}{A character vector of gene names with one value for each row of \code{x}.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{design}{A numeric matrix containing blocking terms, i.e., uninteresting factors driving expression across cells.}
\item{pval.type}{A string specifying the type of combined p-value to be computed, i.e., Simes' or IUT.}
\item{direction}{A string specifying the direction of log-fold changes to be considered for each cluster.}
\item{lfc}{A positive numeric scalar specifying the log-fold change threshold to be tested against.}
\item{log.p}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{full.stats}{A logical scalar indicating whether all statistics (i.e., raw and BH-adjusted p-values) should be returned for each pairwise comparison.}
\item{subset.row}{A logical, integer or character scalar specifying the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to the ANY method.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\details{
This function provides a convenience wrapper for marker gene identification, based on running \code{\link{pairwiseTTests}} and passing the result to \code{\link{combineMarkers}}.
All of the arguments above are supplied directly to one of these two functions.

Note that, if \code{full.stats=TRUE}, the p-values for each pairwise comparison will be log-transformed regardless of the value of \code{log.p}.
}

\value{
A named list of \linkS4class{DataFrame}s, each of which contains a sorted marker gene list for the corresponding cluster.
See \code{?\link{combineMarkers}} for more details on the output format.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{pairwiseTTests}},
\code{\link{combineMarkers}}
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
kout <- kmeans(t(logcounts(y2)), centers=2) # Any clustering method is okay.
out <- findMarkers(y2, clusters=kout$cluster)
}
