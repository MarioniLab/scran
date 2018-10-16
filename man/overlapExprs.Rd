\name{overlapExprs}
\alias{overlapExprs}
\alias{overlapExprs,ANY-method}
\alias{overlapExprs,SingleCellExperiment-method}

\title{Overlap expression profiles}
\description{Compute the gene-specific overlap in expression profiles between two groups of cells.}

\usage{
\S4method{overlapExprs}{ANY}(x, groups, gene.names=rownames(x), block=NULL, 
    pval.type=c("any", "all"), direction=c("any", "up", "down"), tol=1e-8, 
    log.p=FALSE, full.stats=FALSE, subset.row=NULL, BPPARAM=SerialParam()) 

\S4method{overlapExprs}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="logcounts",
    get.spikes=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix of expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{groups}{
A vector of group assignments for all cells.
}
\item{gene.names}{A character vector of gene names with one value for each row of \code{x}.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{pval.type}{A string specifying the type of combined p-value to be computed, i.e., Simes' or IUT.}
\item{direction}{A string specifying which direction of change in expression should be used to rank genes in the output.}
\item{tol}{A numeric scalar specifying the tolerance with which ties are considered.}
\item{log.p}{A logical scalar indicating if log-transformed p-values/FDRs should be returned.}
\item{full.stats}{A logical scalar indicating whether all statistics (i.e., raw and BH-adjusted p-values) should be returned for each pairwise comparison.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{...}{Additional arguments to pass to the matrix method.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{See \code{?"\link{scran-gene-selection}"}.}
}

\details{
This function provides a convenience wrapper for marker gene identification, based on running \code{\link{pairwiseTTests}} and passing the result to \code{\link{combineMarkers}}.
All of the arguments above are supplied directly to one of these two functions.

Note that \code{log.p} only affects the combined p-values and FDRs.
If \code{full.stats=TRUE}, the p-values for each pairwise comparison will be log-transformed regardless of the value of \code{log.p}.
}

\value{
A named list of \linkS4class{DataFrame}s, each of which contains a sorted marker gene list for the corresponding cluster.
See \code{?\link{combineMarkers}} for more details on the output format.
}

\seealso{
See \code{\link{pairwiseWilcox}} and \code{\link{combineMarkers}} for the component functions.

See \code{\link{findMarkers}} for the equivalent function using t-tests.
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
groups <- sample(3, ncol(y2), replace=TRUE)
out <- overlapExprs(y2, groups, subset.row=1:10)
}
