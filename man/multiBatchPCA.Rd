\name{multiBatchPCA}
\alias{multiBatchPCA}

\title{Multi-batch PCA}
\description{Perform a PCA across multiple gene expression matrices to project all cells to a common low-dimensional space.}

\usage{
multiBatchPCA(..., d=50, approximate=FALSE, irlba.args=list(), 
    BPPARAM=SerialParam()) 
}

\arguments{
\item{...}{Two or more matrices containing expression values (usually log-normalized).
Each matrix is assumed to represent one batch.}
\item{d}{An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.}
\item{approximate}{A logical scalar indicating whether \code{\link[irlba]{irlba}} should be used to perform the initial PCA.}
\item{irlba.args}{A list of arguments to pass to \code{\link[irlba]{irlba}} when \code{approximate=TRUE}.}
\item{BPPARAM}{A BiocParallelParam object specifying whether the SVD should be parallelized.}
}

\details{
This function is roughly equivalent to \code{cbind}ing all matrices in \code{...} and performing PCA on the merged matrix.
The difference (aside from greater computational efficiency) is that each sample contributes equally to the identification of the loading vectors.
Specifically, the mean vector used for centering is defined as the grand mean of the mean vectors within each batch.
Each batch's contribution to the gene-gene covariance matrix is also divided by the number of cells.

In effect, we weight the cells in each batch to mimic the situation where all batches have the same number of cells.
This avoids domination of the PCA by batches with a large number of cells, 
and ensures that the variance due to unique subpopulations in smaller batches can be captured.
}

\value{
A list of numeric matrices where each matrix corresponds to a batch and contains the first \code{d} PCs (columns) for all cells in the batch (rows).
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
