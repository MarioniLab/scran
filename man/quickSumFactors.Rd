\name{quickSumFactors}
\alias{quickSumFactors}
\alias{quickSumFactors,ANY-method}
\alias{quickSumFactors,SingleCellExperiment-method}

\title{Quickly compute size factors via summation}
\description{Perform a quick calculation of cell-specific size factors based on summation and normalization against nearest neighbors.}

\usage{
\S4method{quickSumFactors}{ANY}(x, k=20, trend.args=list(), approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, min.mean=1, block=NULL, 
    BNPARAM=NULL, BPPARAM=SerialParam())

\S4method{quickSumFactors}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="counts", 
    get.spikes=FALSE, sf.out=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix-like object of counts, where rows are genes and columns are cells.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{k}{Integer scalar specifying the number of nearest neighbors.}
\item{trend.args}{A list containing named arguments to pass to \code{\link{trendVar}}.}
\item{approximate}{Logical scalar indicating whether an approximate PCA should be performed.}
\item{irlba.args}{List containing arguments to pass to \code{\link[irlba]{irlba}} if \code{approximate=TRUE}.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{min.mean}{A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.}
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
\item Identifies the \code{k} nearest neighbors for each cell, based on the first few PCs from the log-transformed library size-normalized expression matrix.
The number of PCs is identified by running \code{\link{trendVar}} and \code{\link{denoisePCA}}.
\item Computes a library size factor for each cell, normalized against the average pseudo-cell constructed from the neighbor set.
This avoids problems with composition biases, under the assumption that negligible DE is present between neighbors.
\item Rescales the size factor across cells, using a median-based normalization method on the neighbor-derived pseudo-cell for each cell.
The use of a summed profile avoids problems with zeros upon median-based normalization.
}

All cells are rescaled towards a reference pseudo-cell, defined as the cell with the lowest distance to the \code{d}th neighbor.
This is a rough proxy for the neighbourhood with the highest density of cells, which should also have the most stable neighbor-derived pseudo-cell.

If \code{block} is specified, the above process is performed for all cells within each level of block.
To calibrate the size factors across blocks, median-based normalization is applied to the block-specific reference pseudo-cells.
The global reference pseudo-cell is defined as that which is closest to the \dQuote{middle} of the data set.
This is defined as the median position on the first PC after performing a PCA on all of the block-specific reference pseudo-cells.

If \code{min.mean} is positive, genes with library size-normalized averages below \code{min.mean} will be ignored during median-based normalization.
(Specifically, the threshold is applied to the grand average of the pooled profile for each cell and that of the reference cell.)
This favours moderate- and high-abundance genes for which median-based normalization is more accurate.
Note that this filtering will only affect the calculation of scaling factors, and will not be used for nearest-neighbor detection.
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
