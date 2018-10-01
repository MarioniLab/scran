\name{quickSumFactors}
\alias{quickSumFactors}
\alias{quickSumFactors,ANY-method}
\alias{quickSumFactors,SingleCellExperiment-method}

\title{Quickly compute size factors via summation}
\description{Perform a quick calculation of cell-specific size factors based on summation and normalization against nearest neighbors.}

\usage{
\S4method{quickSumFactors}{ANY}(x, k=20, d=50, approximate=FALSE, irlba.args=list(),
    min.mean=1, subset.row=NULL, block=NULL, BNPARAM=NULL, 
    BPPARAM=SerialParam())

\S4method{quickSumFactors}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="counts", 
    get.spikes=FALSE, sf.out=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix-like object of counts, where rows are genes and columns are cells.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{k}{Integer scalar specifying the number of nearest neighbors.}
\item{d}{Integer scalar specifying the number of PCs to retain for nearest-neighbor searching.}
\item{approximate}{Logical scalar indicating whether an approximate PCA should be performed.}
\item{irlba.args}{List containing arguments to pass to \code{\link{irlba}} if \code{approximate=TRUE}.}
\item{min.mean}{A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{block}{An optional factor specifying blocking levels for all cells.}
\item{BNPARAM}{A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for nearest neighbor detection.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying the parallelization to use.}
\item{...}{Additional arguments to pass to \code{quickSumFactors,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{A logical scalar specifying whether spike-in transcripts should be used.}
\item{sf.out}{A logical scalar indicating whether only size factors should be returned.}
}

\details{
This function implements a greatly simplified counterpart to the \code{\link{quickCluster}} and \code{\link{computeSumFactors}} workflow for computing cell-specific size factors.
Briefly, it:
\enumerate{
\item Identifies the \code{k} nearest neighbors for each cell, based on the first \code{d} PCs from the log-transformed library size-normalized expression matrix.
If \code{d} is \code{NA}, the log-expression matrix is used directly without any PCA.
\item Computes a library size factor for each cell, normalized against the average pseudo-cell constructed from the neighbor set.
This avoids problems with composition biases, under the assumption that negligible DE is present between neighbors.
\item Rescales the size factor across cells, using a median-based normalization method on the neighbor-derived pseudo-cell for each cell.
The use of a summed profile avoids problems with zeros upon median-based normalization.
}
}

\value{
For \code{quickSumFactors,ANY-method}, a numeric vector of size factors for all cells in \code{x} is returned.

For \code{quickSumFactors,SingleCellExperiment-method}, an object of class \code{x} is returned containing the vector of size factors in \code{sizeFactors(x)}, if \code{sf.out=FALSE}.
Otherwise, the vector of size factors is returned directly.
}

\author{
Aaron Lun
}

\examples{
# Mocking up some data.
set.seed(100)
popsize <- 200
ngenes <- 10000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs*10, size=1), ncol=popsize, byrow=TRUE)

# Computing the size factors.
out.facs <- quickSumFactors(counts)
head(out.facs)
plot(colSums(counts), out.facs, log="xy")
}

\seealso{
\code{\link{quickCluster}} and \code{\link{computeSumFactors}}, on which this function was based.
}

\keyword{normalization}

\references{
Lun ATL, Bach K and Marioni JC (2016).
Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
\emph{Genome Biol.} 17:75
}
