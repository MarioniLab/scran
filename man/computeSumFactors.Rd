\name{Deconvolution Methods}
\alias{computeSumFactors}
\alias{computeSumFactors,ANY-method}
\alias{computeSumFactors,SingleCellExperiment-method}

\title{Normalization by deconvolution}
\description{Methods to normalize single-cell RNA-seq data by deconvolving size factors from cell pools.}

\usage{
\S4method{computeSumFactors}{ANY}(x, sizes=seq(21, 101, 5), clusters=NULL, ref.clust=NULL,
    max.cluster.size=3000, positive=TRUE, scaling=NULL, min.mean=1,
    subset.row=NULL, BPPARAM=SerialParam())

\S4method{computeSumFactors}{SingleCellExperiment}(x, ..., subset.row=NULL, assay.type="counts",
    get.spikes=FALSE, sf.out=FALSE)
}

\arguments{
\item{x}{
    A numeric matrix-like object of counts, where rows are genes and columns are cells.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{sizes}{A numeric vector of pool sizes, i.e., number of cells per pool.}
\item{clusters}{An optional factor specifying which cells belong to which cluster, for deconvolution within clusters.}
\item{ref.clust}{A level of \code{clusters} to be used as the reference cluster for inter-cluster normalization.}
\item{max.cluster.size}{An integer scalar specifying the maximum number of cells in each cluster.}
\item{positive}{A logical scalar indicating whether linear inverse models should be used to enforce positive estimates.}
\item{scaling}{A numeric scalar containing scaling factors to adjust the counts prior to computing size factors.}
\item{min.mean}{A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{BPPARAM}{A BiocParallelParam object specifying whether and how clusters should be processed in parallel.}
\item{...}{Additional arguments to pass to \code{computeSumFactors,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{See \code{?"\link{scran-gene-selection}"}.}
\item{sf.out}{A logical scalar indicating whether only size factors should be returned.}
}

\section{Overview of the deconvolution method}{
The \code{computeSumFactors} function provides an implementation of the deconvolution strategy (Lun et al., 2016) for scaling normalization of sparse count data.
Briefly, a pool of cells is selected and the expression profiles for those cells are summed together.
The pooled expression profile is normalized against an average reference pseudo-cell, constructed by averaging the counts across all cells.
This defines a size factor for the pool as the median ratio between the count sums and the average across all genes.

The scaling bias for the pool is equal to the sum of the biases for the constituent cells.
The same applies for the size factors, as these are effectively estimates of the bias for each cell.
This means that the size factor for the pool can be written as a linear equation of the size factors for the cells.
Repeating this process for multiple pools will yield a linear system that can be solved to obtain the size factors for the individual cells.

In this manner, pool-based factors are deconvolved to yield the relevant cell-based factors.
The advantage is that the pool-based estimates are more accurate, as summation reduces the number of stochastic zeroes and the associated bias of the size factor estimate.
This accuracy will feed back into the deconvolution process, thus improving the accuracy of the cell-based size factors.
}

\section{Pooling with a sliding window}{
Within each cluster (if not specified, all cells are put into a single cluster), cells are sorted by increasing library size and a sliding window is applied to this ordering.
Each location of the window defines a pool of cells with similar library sizes.
This avoids inflated estimation errors for very small cells when they are pooled with very large cells.
Sliding the window will construct an over-determined linear system that can be solved by least-squares methods to obtain cell-specific size factors.

Window sliding is repeated with different window sizes to construct the linear system, as specified by \code{sizes}.
By default, the number of cells in each window ranges from 21 to 101.
Using a range of window sizes improves the precision of the estimates, at the cost of increased computational complexity.
The defaults were chosen to provide a reasonable compromise between these two considerations.
The default choice also avoids rare cases of linear dependencies and unstable estimates when all pool sizes are not co-prime with the number of cells.

The smallest window should be large enough so that the pool-based size factors are, on average, non-zero.
We recommend window sizes no lower than 20 for UMI data, though smaller windows may be possible for read count data.
}

\section{Prescaling to improve accuracy}{ 
The simplest approach to pooling is to simply add the counts together for all cells in each pool.
However, this is suboptimal as any errors in the estimation of the pooled size factor will propagate to all component cell-specific size factors upon solving the linear system.
If the error is distributed evenly across all cell-specific size factors, the small size factors will have larger relative errors compared to the large size factors.

To avoid this, we perform \dQuote{prescaling} where we divide the counts by a cell-specific factor prior to pooling.
The prescaling factor should be close to the true size factor for each cell.
(Obviously, the true size factor is unknown so we use the library size for each cell instead.)
Solving the linear system constructed with prescaled values should yield estimates that are more-or-less equal across all cells.
Thus, given similar absolute errors, the relative errors for all cells will also be similar.

The default usage assumes that the library size is roughly proportional to the true size factor.
This can be violated in pathological scenarios involving extreme differential expression.
In such cases, we recommend running \code{computeSumFactors} twice to improve accuracy.
The first run is done as usual and will yield an initial estimate of the size factor for each cell.
In the second run, we use our initial estimates for prescaling by supplying them to the \code{scaling} argument, 
This weakens the assumption above by avoiding large inaccuracies in the prescaling factor.

Obviously, this involves twice as much computational work.
As a result, performing multiple iterations is not the default recommendation, especially as the benefits are only apparent in extreme circumstances.
}

\section{Solving the linear system}{
The linear system is solved using the sparse QR decomposition from the \pkg{Matrix} package.
However, this has known problems when the linear system becomes too large (see \url{https://stat.ethz.ch/pipermail/r-help/2011-August/285855.html}).
In such cases, set \code{clusters} to break up the linear system into smaller, more manageable components that can be solved separately.
The default \code{max.cluster.size} will arbitrarily break up the cell population (within each cluster, if specified) so that we never pool more than 3000 cells.
}

\section{Normalization within and between clusters}{
In general, it is more appropriate to pool more similar cells to avoid violating the assumption of a non-DE majority of genes across the data set.
This can be done by specifying the \code{clusters} argument where cells in each cluster have similar expression profiles.
Deconvolution is subsequently applied on the cells within each cluster.
Each cluster should contain a sufficient number of cells for pooling -- a warning is raised if the number of cells is less than the maximum value of \code{sizes}.
A convenince function \code{\link{quickCluster}} is provided for this purpose.

Size factors computed within each cluster must be rescaled for comparison between clusters.
This is done by normalizing between clusters to identify the rescaling factor.
One cluster is chosen as a ``reference'' to which all others are normalized.
Ideally, the reference cluster should have a stable expression profile and not be extremely different from all other clusters.

By default, the cluster with the most non-zero counts is used as the reference.
This improves precision (and also accuracy) of the median-based estimate of the rescaling factor.
Alternatively, the reference can be manually specified using \code{ref.clust} if there is prior knowledge about which cluster is most suitable, e.g., from PCA or t-SNE plots.
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
This results in many zero counts and pooled size factors of zero, and can be fixed by filtering out more genes with a higher \code{min.mean}.
}
Another approach is to increase in the number of \code{sizes} to improve the precision of the estimates.
This reduces the chance of obtaining negative size factors due to estimation error, for cells where the true size factors are very small.

As a last resort, \code{positive=TRUE} is set by default, which uses \code{\link{cleanSizeFactors}} to coerce any negative estimates to positive values.
This ensures that, at the very least, downstream analysis is possible even if the size factors for affected cells are not accurate.
Users can skip this step by setting \code{positive=FALSE} to perform their own diagnostics or coercions.
}

\section{Gene selection}{
Note that pooling does not eliminate the need to filter out low-abundance genes.
As mentioned above, if too many genes have consistently low counts across all cells, even the pool-based size factors will be zero.
This results in zero or negative size factor estimates for many cells.
Filtering ensures that this is not the case, e.g., by removing genes with average counts below 1.

In general, genes with average counts below 1 (for read count data) or 0.1 (for UMI data) should not be used for normalization.
Such genes will automatically be filtered out by applying a minimum threshold \code{min.mean} on the library size-adjusted average counts from \code{\link[scater]{calcAverage}}.
If \code{clusters} is specified, filtering by \code{min.mean} is performed on the per-cluster average during within-cluster normalization,
and then on the (library size-adjusted) average of the per-cluster averages during between-cluster normalization.

Spike-in transcripts are not included with the default \code{get.spikes=FALSE} as they can behave differently from the endogenous genes. 
Users wanting to perform spike-in normalization should see \code{\link{computeSpikeFactors}} instead.
}

\section{Obtaining standard errors}{
Previous versions of \code{computeSumFactors} would return the standard error for each size factor when \code{errors=TRUE}.
This argument is no longer available as we have realized that standard error estimation from the linear model is not reliable.
Errors are likely underestimated due to correlations between pool-based size factors when they are computed from a shared set of underlying counts.
Users wishing to obtain a measure of uncertainty are advised to perform simulations instead, using the original size factor estimates to scale the mean counts for each cell.
Standard errors can then be calculated as the standard deviation of the size factor estimates across simulation iterations.
}

\value{
For \code{computeSumFactors,ANY-method}, a numeric vector of size factors for all cells in \code{x} is returned.

For \code{computeSumFactors,SingleCellExperiment-method}, an object of class \code{x} is returned containing the vector of size factors in \code{sizeFactors(x)}, if \code{sf.out=FALSE}.
Otherwise, the vector of size factors is returned directly.
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
