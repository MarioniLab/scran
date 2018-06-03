\name{fastMNN}
\alias{fastMNN}

\title{Fast mutual nearest neighbors correction}
\description{Correct for batch effects in single-cell expression data using the mutual nearest neighbors method.}

\usage{
fastMNN(..., k=20, cos.norm=TRUE, d=50, ndist=3, approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, BPPARAM=SerialParam()) 
}

\arguments{
\item{...}{Two or more expression matrices where genes correspond to rows and cells correspond to columns.
Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider when identifying mutual nearest neighbors.}
\item{cos.norm}{A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.}
\item{d}{An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.}
\item{ndist}{A numeric scalar specifying the threshold beyond which neighbours are to be ignored when computing correction vectors.
Each threshold is defined in terms of the number of median distances.}
\item{approximate}{A logical scalar indicating whether \code{\link[irlba]{irlba}} should be used to perform the initial PCA.}
\item{irlba.args}{A list of arguments to pass to \code{\link[irlba]{irlba}} when \code{pc.approx=TRUE}.}
\item{subset.row}{A vector specifying the genes to use, providing equivalent behaviour to subsetting each object in \code{...}.
All genes are used by default.}
\item{BPPARAM}{A BiocParallelParam object specifying whether the nearest-neighbor searches should be parallelized.}
}

\value{
A named list containing two components:
\describe{
\item{\code{corrected}:}{A matrix with number of columns equal to \code{d}, and number of rows equal to the total number of cells in \code{...}.
Cells are ordered in the same manner as supplied in \code{...}.}
\item{\code{batch}:}{An integer vector specifying the batch of origin in \code{...} for each row (i.e., cell) in \code{corrected}.}
}}

\details{
This function provides a variant of the \code{\link{mnnCorrect}} function, modified for speed and more robust performance.
In particular, it performs a multi-sample PCA prior to any correction, and performs all calculations in the PC space.
This reduces computational work and provides some denoising (see \code{\link{denoisePCA}} for details).
As a result, though, the corrected output cannot be interpreted on a gene level and is useful only for cell-level comparisons, e.g., clustering and visualization.

Another simplification is that correction vector for each cell is directly computed from its \code{k} nearest neighbours in the same batch.
Specifically, only the \code{k} nearest neighbouring cells that \emph{also} participate in MNN pairs are used.
Each MNN-participating neighbour is weighted by distance from the current cell, using a tricube scheme with bandwidth equal to the median distance multiplied by \code{ndist}.
This ensures that the correction vector only uses information from the closest cells, improving the fidelity of local correction.

Finally, issues with \dQuote{kissing} are avoided with a two-step procedure that removes variation along the batch effect vector.
First, the average correction vector across all MNN pairs is computed.
Cell coordinates are adjusted such that all cells in a single batch have the same position along this vector.
The correction vectors are then recalculated with the adjusted coordinates (but the same MNN pairs).
}

\author{
Aaron Lun    
}

\seealso{
\code{\link{mnnCorrect}},
\code{\link[irlba]{irlba}}
}

\examples{
B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
out <- fastMNN(B1, B2) # corrected values
}
