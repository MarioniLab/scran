\name{simpleSumFactors}
\alias{simpleSumFactors}
\alias{simpleSumFactors,ANY-method}
\alias{simpleSumFactors,SingleCellExperiment-method}

\title{Quickly compute size factors via summation}
\description{Perform a simple calculation of cell-specific size factors based on summation and normalization against nearest neighbors.}

\usage{
\S4method{simpleSumFactors}{ANY}(x, k=20, ref.k=k, trend.args=list(), approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, min.mean=1, block=NULL, 
    BNPARAM=NULL, BPPARAM=SerialParam())

\S4method{simpleSumFactors}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="counts", 
    get.spikes=FALSE, sf.out=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix-like object of counts, where rows are genes and columns are cells.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{k}{Integer scalar specifying the number of nearest neighbors.}
\item{ref.k}{Integer scalar specifying the number of nearest neighbors to use for constructing the reference pseudo-cell.}
\item{trend.args}{A list containing named arguments to pass to \code{\link{trendVar}}.}
\item{approximate}{Logical scalar indicating whether an approximate PCA should be performed.}
\item{irlba.args}{List containing arguments to pass to \code{\link[irlba]{irlba}} if \code{approximate=TRUE}.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{min.mean}{A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.}
\item{block}{An optional factor specifying blocking levels for all cells.}
\item{ref.block}{An integer or string specifying the level of \code{block} to use as the reference.}
\item{BNPARAM}{A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for nearest neighbor detection.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying the parallelization to use.}
\item{...}{Additional arguments to pass to \code{simpleSumFactors,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{See \code{?"\link{scran-gene-selection}"}.}
\item{sf.out}{A logical scalar indicating whether only size factors should be returned.}
}

\details{
This function implements a greatly simplified counterpart to the \code{\link{quickCluster}} and \code{\link{computeSumFactors}} workflow for computing cell-specific size factors.
Briefly, it:
\enumerate{
\item Identifies the \code{k} nearest neighbors for each cell, based on the first few PCs from the log-transformed library size-normalized expression matrix.
The number of PCs is determined running \code{\link{trendVar}} and \code{\link{denoisePCA}}.
\item Constructs a pseudo-cell from the neighbor set by taking the weighted average of the neighbors' expression profiles.
Weights are computed using a tricube kernel on the distance, with bandwidth equal to three-fold the median distance.
This protects against large distances for cell types with fewer than \code{k} cells in the data set.
\item Computes the ratio of the library size for each cell to the sum of averaged values for the associated pseudo-cell.
This represents an initial size factor estimate, equivalent to library size normalization.
The assumption is that neighboring cells are from the same cell type and do not exhibit composition biases associated with DE.
\item Performs median-based normalization on the neighbor-derived pseudo-cell for each cell to a reference pseudo-cell.
This yields a rescaling factor that is multiplied by the initial ratio to obtain the cell-specific size factor.
The use of the median protects against DE between pseudo-cells, while the use of a summed profile avoids problems with zeros.
}

The number of neighbors \code{k} represents the minimum number of cells of any cell type or state.
Reducing \code{k} increases the robustness of the algorithm to population heterogeneity.
However, it reduces the precision of the median-based estimator and increases the risk of obtaining zero or undefined size factors after rescaling.
\code{k} should be set to the largest tolerable value for optimal performance, though the default is usually sufficient.

The reference pseudo-cell is constructed from the neighbors of the cell with the lowest distance to the \code{ref.k}th neighbor.
This is a rough proxy for the neighbourhood with the highest density of cells, which should also have the most stable neighbor-derived pseudo-cell.
By default, \code{ref.k} is equal to \code{k} but can be increased to provide a more stable reference pseudo-cell.
(This is more specific than increasing \code{k}, which would require stronger assumptions for \emph{all} cells.)

If \code{block} is specified, the above process is performed separately for all cells within each level of \code{block}.
To calibrate the size factors across blocks, median-based normalization is applied to the block-specific reference pseudo-cells to obtain a calibration factor.
The global reference pseudo-cell is specified with \code{ref.block} or, if \code{NULL}, is defined as that which is closest to the \dQuote{middle} of the data set.
This refers to the median position on the first PC after performing a PCA on the log-transformed expression profiles of the block-specific reference pseudo-cells.
The size factor for each cell in each block is then multiplied by the calibration factor for that block to obtain the \emph{final} size factor for each cell. 

If \code{min.mean} is positive, genes with library size-normalized averages below \code{min.mean} will be ignored during median-based normalization.
Specifically, the threshold is applied to the grand average of the pooled profile for each cell and that of the reference cell.
This favours moderate- and high-abundance genes for which median-based normalization is more accurate.
Note that this filtering will only affect the calculation of scaling and calibration factors, and will not be used for nearest-neighbor detection.
}

\section{Trouble-shooting errors and warnings}{
\dQuote{cleaning zero size factor estimates}: the rescaling factor is zero for some pseudo-cells. 
Increase \code{k} to reduce the number of stochastic zeroes in each pseudo-cell.
 
\dQuote{infinite size factors}: the rescaling factor is \code{Inf} for some pseudo-cells.
Increase \code{ref.k} to reduce the number of stochastic zeroes in the reference pseudo-cell.

\dQuote{no genes remaining after filtering}: filtering has removed all genes for some pseudo-cells.
Lower \code{min.mean} to retain more genes.

\dQuote{calibration factors are not strictly positive}: some of the inter-block calibration factors are not positive.
Increase \code{ref.k} to reduce the number of stochastic zeroes in each block-specific reference psuedo-cell.
}

\value{
For \code{simpleSumFactors,ANY-method}, a numeric vector of size factors for all cells in \code{x} is returned.

For \code{simpleSumFactors,SingleCellExperiment-method}, an object of class \code{x} is returned containing the vector of size factors in \code{sizeFactors(x)}, if \code{sf.out=FALSE}.
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
out.facs <- simpleSumFactors(counts)
head(out.facs)
plot(colSums(counts), out.facs, log="xy")
}

\seealso{
\code{\link{quickCluster}} and \code{\link{computeSumFactors}}, on which this function was based.
}

\keyword{normalization}

\references{
Lun ATL (2018).
Further scaling normalization.
\url{https://github.com/MarioniLab/FurtherNorm2018}

Lun ATL, Bach K and Marioni JC (2016).
Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
\emph{Genome Biol.} 17:75
}
