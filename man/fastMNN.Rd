\name{fastMNN}
\alias{fastMNN}

\title{Fast mutual nearest neighbors correction}
\description{Correct for batch effects in single-cell expression data using the mutual nearest neighbors (MNN) method.}

\usage{
fastMNN(..., k=20, cos.norm=TRUE, d=50, ndist=3, approximate=FALSE, 
    irlba.args=list(), auto.order=FALSE,
    BPPARAM=SerialParam()) 
}

\arguments{
\item{...}{Two or more log-expression matrices where genes correspond to rows and cells correspond to columns.
Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.}
\item{cos.norm}{A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.}
\item{d}{An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.}
\item{ndist}{A numeric scalar specifying the threshold beyond which neighbours are to be ignored when computing correction vectors.
Each threshold is defined in terms of the number of median distances.}
\item{approximate}{A logical scalar indicating whether \code{\link[irlba]{irlba}} should be used to perform the initial PCA.}
\item{irlba.args}{A list of arguments to pass to \code{\link[irlba]{irlba}} when \code{pc.approx=TRUE}.}
\item{auto.order}{Logical scalar indicating whether re-ordering of batches should be performed to maximize the number of MNN pairs at each step.}
\item{BPPARAM}{A BiocParallelParam object specifying whether the nearest-neighbor searches should be parallelized.}
}

\value{
A named list containing:
\describe{
\item{\code{corrected}:}{A matrix with number of columns equal to \code{d}, and number of rows equal to the total number of cells in \code{...}.
Cells are ordered in the same manner as supplied in \code{...} if \code{auto.order=FALSE},
otherwise the batches will be re-ordered as determined by the algorithm.}
\item{\code{origin}:}{A DataFrame containing the \code{batch} field (integer, or character if \code{...} is named), 
which specifies the batch of origin for each row (i.e., cell) in \code{corrected};
and the \code{cell} field (integer), which specifies the column of origin for each row in \code{corrected}.} 
\item{\code{pairs}:}{A list of DataFrames specifying which pairs of cells in \code{corrected} were identified as MNNs at each step.} 
}}

\details{
This function provides a variant of the \code{\link{mnnCorrect}} function, modified for speed and more robust performance.
In particular:
\itemize{
\item It performs a multi-sample PCA prior to any correction, and performs all calculations in the PC space.
This reduces computational work and provides some denoising (see \code{\link{denoisePCA}} for details).
As a result, though, the corrected output cannot be interpreted on a gene level and is useful only for cell-level comparisons, e.g., clustering and visualization.
\item The correction vector for each cell is directly computed from its \code{k} nearest neighbours in the same batch.
Specifically, only the \code{k} nearest neighbouring cells that \emph{also} participate in MNN pairs are used.
Each MNN-participating neighbour is weighted by distance from the current cell, using a tricube scheme with bandwidth equal to the median distance multiplied by \code{ndist}.
This ensures that the correction vector only uses information from the closest cells, improving the fidelity of local correction.
\item Issues with \dQuote{kissing} are avoided with a two-step procedure that removes variation along the batch effect vector.
First, the average correction vector across all MNN pairs is computed.
Cell coordinates are adjusted such that all cells in a single batch have the same position along this vector.
The correction vectors are then recalculated with the adjusted coordinates (but the same MNN pairs).
}

By default, batches are merged in the user-supplied order.
However, if \code{auto.order=TRUE}, batches are ordered to maximize the number of MNN pairs at each step.
The aim is to improve the stability of the correction by first merging more similar batches with more MNN pairs.
This can be somewhat time-consuming as MNN pairs need to be iteratively recomputed for all possible batch pairings,
and it is often more convenient for the user to specify an appropriate ordering based on prior knowledge about the batches.

The default setting of \code{cos.norm=TRUE} provides some protection against differences in scaling for arbitrary log-expression matrices.
If possible, we recommend equalizing coverage at the count level before log-transformation - see \code{\link{multiBatchNorm}} for details. 
}

\author{
Aaron Lun    
}

\seealso{
\code{\link{mnnCorrect}},
\code{\link[irlba]{irlba}},
\code{\link{multiBatchNorm}}
}

\examples{
B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
out <- fastMNN(B1, B2) # corrected values
}
