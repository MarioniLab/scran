\name{Deconvolution Methods}
\alias{computeSumFactors}
\alias{computeSumFactors,matrix-method}
\alias{computeSumFactors,SCESet-method}

\title{Normalization by deconvolution}
\description{Methods to normalize single-cell RNA-seq data by deconvolving size factors from cell pools.}

\usage{
\S4method{computeSumFactors}{matrix}(x, sizes=seq(20, 100, 5), clusters=NULL, 
    ref.clust=NULL, positive=FALSE, errors=FALSE, subset.row=NULL)
\S4method{computeSumFactors}{SCESet}(x, subset.row=NULL, ..., assay="counts", 
    get.spikes=FALSE, sf.out=FALSE)
}

\arguments{
\item{x}{
    A numeric count matrix where rows are genes and columns are cells.
    Alternatively, a SCESet object containing such a matrix.
}
\item{sizes}{A numeric vector of pool sizes, i.e., number of cells per pool.}
\item{clusters}{An optional factor specifying which cells belong to which cluster, for deconvolution within clusters.}
\item{ref.clust}{A level of \code{clusters} to be used as the reference cluster for inter-cluster normalization.}
\item{positive}{A logical scalar indicating whether linear inverse models should be used to enforce positive estimates.}
\item{errors}{A logical scalar indicating whether the standard error should be returned.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to \code{computeSumFactors,matrix-method}.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical scalar specifying whether spike-in transcripts should be used.}
\item{sf.out}{A logical scalar indicating whether only size factors should be returned.}
}

\section{Overview of the deconvolution method}{
The \code{computeSumFactors} function provides an implementation of the deconvolution strategy for normalization.
Briefly, a pool of cells is selected and the counts for those cells are summed together.
The count sums for this pool is normalized against an average reference pseudo-cell, constructed by averaging the counts across all cells.
This defines a size factor for the pool as the median ratio between the count sums and the average across all genes.

Now, the bias for the pool is equal to the sum of the biases for the constituent cells.
The same applies for the size factors (which are effectively estimates of the bias for each cell).
This means that the size factor for the pool can be written as a linear equation of the size factors for the cells.
Repeating this process for multiple pools will yield a linear system that can be solved to obtain the size factors for the individual cells.

In this manner, pool-based factors are deconvolved to yield the relevant cell-based factors.
The advantage is that the pool-based estimates are more accurate, as summation reduces the number of stochastic zeroes and the associated bias of the size factor estimate.
This accuracy will feed back into the deconvolution process, thus improving the accuracy of the cell-based size factors.
}

\section{Normalization within and between clusters}{
In general, it is more appropriate to pool more similar cells to avoid violating the assumption of a non-DE majority of genes across the data set.
This can be done by specifying the \code{clusters} argument where cells in each cluster have similar expression profiles.
Deconvolution is subsequently applied on the cells within each cluster.
Each cluster should contain a sufficient number of cells for pooling -- an error is thrown if the number of cells is less than the maximum value of \code{sizes}.
A convenince function \code{\link{quickCluster}} is provided for rapid clustering based on Spearman's rank correlation.

Size factors computed within each cluster must be rescaled for comparison between clusters.
This is done by normalizing between clusters to identify the rescaling factor.
One cluster is chosen as a ``reference'' (by default, that with the median of the mean per-cell library sizes is used) to which all others are normalized.
Ideally, a cluster that is not extremely different from all other clusters should be used as the reference.
This can be specified using \code{ref.clust} if there is prior knowledge about which cluster is most suitable, e.g., from PCA or t-SNE plots.
}

\section{Additional details about pooling and deconvolution}{
Within each cluster (if not specified, all cells are put into a single cluster), cells are sorted by increasing library size and a sliding window is applied to this ordering.
Each location of the window defines a pool of cells with similar library sizes.
This avoids inflated estimation errors for very small cells when they are pooled with very large cells.
Sliding the window will construct an over-determined linear system that can be solved with methods like the QR decomposition.

Window sliding is repeated with different window sizes to construct the linear system, as specified by \code{sizes}.
By default, the number of cells in each window ranges from 20 to 100.
Using a range of window sizes improves the precision of the estimates, at the cost of increased computational complexity.
The defaults are chosen to provide a reasonable compromise between these two considerations.

The smallest window should be large enough so that the pool-based size factors are, on average, non-zero.
We recommend window sizes no lower than 20 for UMI data, though smaller windows may be possible for read count data.
}

\section{Dealing with negative size factors}{ 
In theory, it is possible to obtain negative estimates for the size factors.
These values are obviously nonsensical and \code{computeSumFactors} will raise a warning if they are encountered.
Negative estimates are mostly commonly generated from low quality cells with few expressed features, such that most counts are zero even after pooling.
They may also occur if insufficient filtering of low-abundance genes was performed.

To avoid negative size factors, the best solution is to increase the stringency of the filtering.
\itemize{
\item If only a few negative size factors are present, they are likely to correspond to a few low-quality cells with few expressed features.
Such cells are difficult to normalize reliably under any approach, and can be removed by increasing the stringency of the quality control.
\item If many negative size factors are present, it is probably due to insufficient filtering of low-abundance genes.
This results in many zero counts and pooled size factors of zero, and can be fixed by filtering out more genes.
}
A complementary approach is to increase in the number of \code{sizes} to improve the precision of the estimates.
This reduces the chance of obtaining negative size factors due to estimation error for cells with near-zero true size factors.

As a last resort, some protection can be provided by setting \code{positive=TRUE}, which will use linear inverse models to solve the system.
This ensures that non-negative values for the size factors will always be obtained.
Note that some cells may still have size factors of zero and should be removed prior to downstream analysis.
Such occurrences are unavoidable -- rather, the aim is to prevent negative values from affecting the estimates for all other cells.
}

\section{Gene selection}{
By default, \code{get.spikes=FALSE} in \code{quickCluster,SCESet-method} which means that spike-in transcripts are not included in the set of genes used for deconvolution.
This is because they can behave differently from the endogenous genes. 
Users wanting to perform spike-in normalization should see \code{\link{computeSpikeFactors}} instead.

Users can also set \code{subset.row} to specify which rows of \code{x} are to be used to calculate correlations.
This is equivalent to but more efficient than subsetting \code{x} directly, as it avoids constructing a (potentially large) temporary matrix.
If this is specified, it will overwrite any setting of \code{get.spikes}.

Note that pooling does not eliminate the need to filter out low-abundance genes.
As mentioned above, if too many genes have consistently low counts across all cells, even the pool-based size factors will be zero.
This results in zero or negative size factor estimates for many cells.
Filtering ensures that this is not the case, e.g., by removing genes with average counts below 1.
}

\value{
For \code{computeSumFactors,matrix-method}, a numeric vector of size factors for all cells in \code{x} is returned.

For \code{computeSumFactors,SCESet-method}, an object of class \code{x} is returned containing the vector of size factors in \code{sizeFactors(x)}, if \code{sf.out=FALSE}.
Otherwise, the vector of size factors is returned directly.

If \code{errors=TRUE}, the standard errors of the size factor estimates are stored as the \code{"standard.error"} field of the attributes of the returned vector.
}

\author{
Aaron Lun and Karsten Bach
}

\seealso{
\code{\link{quickCluster}}
}

\examples{
# Mocking up some data.
set.seed(100)
popsize <- 200
ngenes <- 10000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs*10, size=1), ncol=popsize, byrow=TRUE)

# Computing the size factors.
out.facs <- computeSumFactors(counts)
head(out.facs)
plot(colSums(counts), out.facs, log="xy")
}

\keyword{normalization}

\references{
Lun ATL, Bach K and Marioni JC (2016).
Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
\emph{Genome Biol.} 17:75
}
