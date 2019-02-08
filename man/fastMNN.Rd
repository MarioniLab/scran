\name{fastMNN}
\alias{fastMNN}

\title{Fast mutual nearest neighbors correction}
\description{Correct for batch effects in single-cell expression data using the mutual nearest neighbors (MNN) method.}

\usage{
fastMNN(..., k=20, cos.norm=TRUE, ndist=3, d=50, approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, auto.order=FALSE, pc.input=FALSE,
    compute.variances=FALSE, assay.type="logcounts", get.spikes=FALSE, 
    BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
}

\arguments{
\item{...}{Two or more log-expression matrices where genes correspond to rows and cells correspond to columns.
One matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).

Alternatively, two or more \linkS4class{SingleCellExperiment} objects can be supplied, where each object contains a log-expression matrix in the \code{assay.type} assay.

Alternatively, two or more matrices of low-dimensional representations can be supplied if \code{pc.input=TRUE}.
Here, rows are cells and columns are dimensions (the latter should be common across all batches).}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.}
\item{cos.norm}{A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.}
\item{ndist}{A numeric scalar specifying the threshold beyond which neighbours are to be ignored when computing correction vectors.
Each threshold is defined in terms of the number of median distances.}
\item{d, approximate, irlba.args}{Further arguments to pass to \code{\link{multiBatchPCA}}.
Setting \code{approximate=TRUE} is recommended for large data sets with many cells.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{auto.order}{Logical scalar indicating whether re-ordering of batches should be performed to maximize the number of MNN pairs at each step.

Alternatively, an integer vector containing a permutation of \code{1:N} where \code{N} is the number of batches.}
\item{compute.variances}{Logical scalar indicating whether the percentage of variance lost due to non-orthogonality should be computed.}
\item{pc.input}{Logical scalar indicating whether the values in \code{...} are already low-dimensional, e.g., the output of \code{\link{multiBatchPCA}}.}
\item{assay.type}{A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.}
\item{get.spikes}{See \code{?"\link{scran-gene-selection}"}. Only relevant if \code{...} contains SingleCellExperiment objects.}
\item{BNPARAM}{A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm to use.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying whether the PCA and nearest-neighbor searches should be parallelized.}
}

\value{
A named list containing:
\describe{
\item{\code{corrected}:}{A matrix with number of columns equal to \code{d}, and number of rows equal to the total number of cells in \code{...}.
Cells are ordered in the same manner as supplied in \code{...}.}
\item{\code{batch}:}{A \linkS4class{Rle} containing the batch of origin for each row (i.e., cell) in \code{corrected}.}
\item{\code{pairs}:}{A list of DataFrames specifying which pairs of cells in \code{corrected} were identified as MNNs at each step.} 
}

If \code{pc.input=FALSE}, an additional \code{rotation} field will be present, containing the matrix of rotation vectors from \code{\link{multiBatchPCA}}.

If \code{compute.variances=TRUE}, an additional \code{lost.var} field will be present containing a numeric vector of percentages of variances lost for each batch.
}

\details{
This function provides a variant of the \code{\link{mnnCorrect}} function, modified for speed and more robust performance.
In particular:
\itemize{
\item It performs a multi-sample PCA via \code{\link{multiBatchPCA}} and subsequently performs all calculations in the PC space.
This reduces computational work and provides some denoising - see, comments in \code{?\link{denoisePCA}}.
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

The default setting of \code{cos.norm=TRUE} provides some protection against differences in scaling for arbitrary expression matrices.
However, if possible, we recommend using the output of \code{\link{multiBatchNorm}} as input to \code{fastMNN}.
This will equalize coverage on the count level before the log-transformation, which is a more accurate rescaling than cosine normalization on the log-values.

If \code{compute.variances=TRUE}, the function will compute the percentage of variance that is parallel to the average correction vectors at each merge step.
This represents the variance that is not orthogonal to the batch effect and subsequently lost upon correction.
Large proportions suggest that there is biological structure that is parallel to the batch effect, 
corresponding to violations of the assumption that the batch effect is orthogonal to the biological subspace.
}

\section{Controlling the merge order}{
By default, batches are merged in the user-supplied order.
However, if \code{auto.order=TRUE}, batches are ordered to maximize the number of MNN pairs at each step.
The aim is to improve the stability of the correction by first merging more similar batches with more MNN pairs.
This can be somewhat time-consuming as MNN pairs need to be iteratively recomputed for all possible batch pairings.
It is often more convenient for the user to specify an appropriate ordering based on prior knowledge about the batches.

If \code{auto.order} is an integer vector, it is treated as an ordering permutation with which to merge batches.
For example, if \code{auto.order=c(4,1,3,2)}, batches 4 and 1 in \code{...} are merged first, followed by batch 3 and then batch 2.
This is often more convenient than changing the order manually in \code{...}, which would alter the order of batches in the output \code{corrected} matrix.
Indeed, no matter what the setting of \code{auto.order} is, the order of cells in the output corrected matrix is \emph{always} the same.

Further control of the merge order can be achieved by performing the multi-sample PCA outside of this function with \code{\link{multiBatchPCA}}.
Then, batches can be progressively merged by repeated calls to \code{fastMNN} with \code{pc.input=TRUE}.
This is useful in situations where the order of batches to merge is not straightforward, e.g., involving hierarchies of batch similarities. 
We only recommend this mode for advanced users, and note that:
\itemize{
    \item \code{\link{multiBatchPCA}} will not perform cosine-normalization, 
so it is the responsibility of the user to cosine-normalize each batch beforehand with \code{\link{cosineNorm}} to recapitulate results with \code{cos.norm=TRUE}.
    \item \code{\link{multiBatchPCA}} must be run on all samples at once, to ensure that all cells are projected to the same low-dimensional space.
    \item Setting \code{pc.input=TRUE} is criticial to avoid unnecessary (and incorrect) cosine-normalization and PCA within each step of the merge.
}
See the Examples below for how the \code{pc.input} argument should be used.
}

\section{Choice of genes}{
Users should set \code{subset.row} to subset the inputs to highly variable genes or marker genes.
This provides more meaningful identification of MNN pairs by reducing the noise from irrelevant genes.
Note that users should not be too restrictive with subsetting, as high dimensionality is required to satisfy the orthogonality assumption in MNN detection.

For SingleCellExperiment inputs, spike-in transcripts should be the same across all objects.
They are automatically removed unless \code{get.spikes=TRUE}, see \code{?"\link{scran-gene-selection}"} for more details.

The reported coordinates for cells refer to a low-dimensional space, but it may be desirable to obtain corrected gene expression values, e.g., for visualization.
This can be done by computing the cross-product of the coordinates with the rotation matrix - see the Examples below.
Note that this will represent corrected values in the space defined by the inputs (e.g., log-transformed) and after cosine normalization if \code{cos.norm=TRUE}.
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
names(out)

# An equivalent approach with PC input.
cB1 <- cosineNorm(B1)
cB2 <- cosineNorm(B2)
pcs <- multiBatchPCA(cB1, cB2)
out.2 <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
all.equal(out[1:3], out.2) # should be TRUE (no rotation)

# Obtaining corrected expression values for genes 1 and 10.
cor.exp <- tcrossprod(out$rotation[c(1,10),], out$corrected)
dim(cor.exp)
}

\references{
Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018).
Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.
\emph{Nat. Biotechnol.} 36(5):421

Lun ATL (2018).
Further MNN algorithm development.
\url{https://github.com/MarioniLab/FurtherMNN2018}
}
