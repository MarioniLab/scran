\name{multiBatchPCA}
\alias{multiBatchPCA}

\title{Multi-batch PCA}
\description{Perform a PCA across multiple gene expression matrices to project all cells to a common low-dimensional space.}

\usage{
multiBatchPCA(..., d=50, approximate=FALSE, irlba.args=list(), 
    subset.row=NULL, assay.type="logcounts", use.spikes=FALSE, 
    BPPARAM=SerialParam()) 
}

\arguments{
\item{...}{Two or more matrices containing expression values (usually log-normalized).
Each matrix is assumed to represent one batch.
Alternatively, two or more SingleCellExperiment objects containing these matrices.}
\item{d}{An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.}
\item{approximate}{A logical scalar indicating whether \code{\link[irlba]{irlba}} should be used to perform the initial PCA.}
\item{irlba.args}{A list of arguments to pass to \code{\link[irlba]{irlba}} when \code{approximate=TRUE}.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{assay.type}{A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.}
\item{use.spikes}{A logical scalar specifying whether spike-in transcripts should be used, if SingleCellExperiment objects are present in \code{...}.}
\item{BPPARAM}{A BiocParallelParam object specifying whether the SVD should be parallelized.}
}

\details{
This function is roughly equivalent to \code{cbind}ing all matrices in \code{...} and performing PCA on the merged matrix.
The difference (aside from greater computational efficiency) is that each sample contributes equally to the identification of the loading vectors.
Specifically, the mean vector used for centering is defined as the grand mean of the mean vectors within each batch.
Each batch's contribution to the gene-gene covariance matrix is also divided by the number of cells.

In effect, we weight the cells in each batch to mimic the situation where all batches have the same number of cells.
This avoids ensures that the variance due to unique subpopulations in smaller batches can be captured.
Otherwise, batches with a large number of cells would dominate the PCA; the mean vector and covariance matrix would be almost fully defined by those batches.

If \code{...} contains SingleCellExperiment objects, any spike-in transcripts should be the same across all batches.
These will be removed prior to PCA unless \code{use.spikes=TRUE}.
Users can also force the algorithm to only use a subset of genes with \code{subset.row}, which will intersect with the non-spike-in genes unless \code{use.spikes=TRUE}.
}

\value{
A \linkS4class{List} of numeric matrices where each matrix corresponds to a batch and contains the first \code{d} PCs (columns) for all cells in the batch (rows).

The metadata also contains a matrix of rotation vectors, which can be used to construct a low-rank approximation of the input matrices.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{fastMNN}}
}

\examples{
d1 <- matrix(rnorm(5000), ncol=100)
d1[1:10,1:10] <- d1[1:10,1:10] + 2 # unique population in d1
d2 <- matrix(rnorm(2000), ncol=40)
d2[11:20,1:10] <- d2[11:20,1:10] + 2 # unique population in d2

out <- multiBatchPCA(d1, d2)

xlim <- range(c(out[[1]][,1], out[[2]][,1]))
ylim <- range(c(out[[1]][,2], out[[2]][,2]))
plot(out[[1]][,1], out[[1]][,2], col="red", xlim=xlim, ylim=ylim)
points(out[[2]][,1], out[[2]][,2], col="blue") 
}
