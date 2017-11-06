\name{mnnCorrect}
\alias{mnnCorrect}

\title{Mutual nearest neighbors correction}
\description{Correct for batch effects in single-cell expression data using the mutual nearest neighbors method.}

\usage{
mnnCorrect(..., k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE,
    svd.dim=0L, subset.row=NULL,  order=NULL, pc.approx=FALSE, 
    BPPARAM=SerialParam())
}

\arguments{
\item{...}{Two or more expression matrices where genes correspond to rows and cells correspond to columns.
Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider when identifying mutual nearest neighbors.}
\item{sigma}{A numeric scalar specifying the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell.}
\item{cos.norm.in}{A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.}
\item{cos.norm.out}{A logical scalar indicating whether cosine normalization should be performed prior to computing corrected expression values.}
\item{svd.dim}{An integer scalar specifying the number of dimensions to use for summarizing biological substructure within each batch.}
\item{subset.row}{A vector specifying the genes with which distances between cells are calculated, e.g., for identifying mutual nearest neighbours.
All genes are used by default.}
\item{order}{An integer vector specifying the order in which batches are to be corrected.}
\item{pc.approx}{A logical scalar indicating whether \code{\link[irlba]{irlba}} should be used to identify the biological subspace.}
\item{BPPARAM}{A BiocParallelParam object specifying whether the nearest-neighbor searches should be parallelized.}
}

\value{
A named list containing two components:
\describe{
\item{\code{corrected}:}{A list of length equal to the number of batches, containing matrices of corrected expression values for each cell in each batch.
The order of batches is the same as supplied in \code{...}, and the order of cells in each matrix is also unchanged.
}
\item{\code{mnn}:}{An integer matrix with two columns, specifying the number of MNNs used for correction of each batch.
Each row corresponds to a batch, and the first and second columns contain the number of MNNs in the reference and current batch, respectively.
}
}}

\details{
This function is designed for batch correction of single-cell RNA-seq data where the batches are partially confounded with biological conditions of interest.
It does so by identifying pairs of mutual nearest neighbors (MNN) in the high-dimensional expression space.
Each MNN pair represents cells in different batches that are of the same cell type/state, assuming that batch effects are mostly orthogonal to the biological manifold.
Correction vectors are calculated from the pairs of MNNs and corrected expression values are returned for use in clustering and dimensionality reduction.

The concept of a MNN pair can be explained by considering cells in each of two batches.
For each cell in one batch, the set of \code{k} nearest cells in the other batch is identified, based on the Euclidean distance in expression space.
Two cells in different batches are considered to be MNNs if each cell is in the other's set.
The size of \code{k} can be interpreted as the minimum size of a subpopulation in each batch.
The algorithm is generally robust to the choice of \code{k}, though values that are too small will not yield enough MNN pairs, while values that are too large will ignore substructure within each batch.

For each MNN pair, a pairwise correction vector is computed based on the difference in the expression profiles.
The correction vector for each cell is computed by applying a Gaussian smoothing kernel with bandwidth \code{sigma} is the pairwise vectors.
This stabilizes the vectors across many MNN pairs and extends the correction to those cells that do not have MNNs.
The choice of \code{sigma} determines the extent of smoothing - a value of 1 is used by default to reflect the boundaries of the space after cosine normalization.
}

\section{Choosing the gene set}{
Distances between cells are calculated with all genes if \code{subset.row=NULL}.
However, users can set \code{subset.row} to perform the distance calculation on a subset of genes, e.g., highly variable genes or marker genes.
This may provide more meaningful identification of MNN pairs by reducing the noise from irrelevant genes.

Regardless of whether \code{subset.row} is specified, corrected values are returned for \emph{all} genes.
This is possible as \code{subset.row} is only used to identify the MNN pairs and other cell-based distance calculations.
Correction vectors between MNN pairs can then be computed in the original space involving all genes in the supplied matrices.
}

\section{Expected type of input data}{
The input expression values should generally be log-transformed, e.g., log-counts, see \code{\link[scater]{normalize}} for details.
They should also be normalized within each data set to remove cell-specific biases in capture efficiency and sequencing depth.
By default, a further cosine normalization step is performed on the supplied expression data to eliminate gross scaling differences between data sets.
\itemize{
\item When \code{cos.norm.in=TRUE}, cosine normalization is performed on the matrix of expression values used to compute distances between cells.
This can be turned off when there are no scaling differences between data sets. 
\item When \code{cos.norm.out=TRUE}, cosine normalization is performed on the matrix of values used to calculate correction vectors (and on which those vectors are applied).
This can be turned off to obtain corrected values on the log-scale, similar to the input data.
}
}

\section{Further options}{
The function depends on a shared biological manifold, i.e., one or more cell types/states being present in multiple batches.
If this is not true, MNNs may be incorrectly identified, resulting in over-correction and removal of interesting biology.
Some protection can be provided by removing components of the correction vectors that are parallel to the biological subspaces in each batch.
The biological subspace in each batch is identified with a SVD on the expression matrix, using either \code{\link{svd}} or \code{\link[irlba]{irlba}}.
The number of dimensions of this subspace can be controlled with \code{svd.dim}.
(By default, this option is turned off by setting \code{svd.dim=0}.)

Users should note that the order in which batches are corrected will affect the final results.
The first batch in \code{order} is used as the reference batch against which the second batch is corrected.
Corrected values of the second batch are added to the reference batch, against which the third batch is corrected, and so on.
This strategy maximizes the chance of detecting sufficient MNN pairs for stable calculation of correction vectors in subsequent batches.
We would consider 20 cells involved in MNN pairs to be the minimum number required for batch correction.
}

\author{
Laleh Haghverdi,
with modifications by Aaron Lun
}

\seealso{
\code{\link[FNN]{get.knnx}}
\code{\link[irlba]{irlba}}
}

\examples{
B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
out <- mnnCorrect(B1, B2) # corrected values
}
