\name{mnnCorrect}
\alias{mnnCorrect}

\title{Mutual nearest neighbors correction}
\description{Correct for batch effects in single-cell expression data using the mutual nearest neighbors method.}

\usage{
mnnCorrect(..., inquiry.genes=NULL, hvg.genes=NULL, k=20, sigma=0.1, 
    cos.norm=TRUE, svd.dim=2, order=NULL) 
}

\arguments{
\item{...}{Two or more expression matrices where genes correspond to rows and cells correspond to columns.
Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).}
\item{inquiry.genes}{A vector specifying the genes in all batches on which the correction is to be performed. 
All genes are used by default.}
\item{hvg.genes}{A vector specifying the genes in all batches with which to identify the mutual nearest neighbours. 
All genes are used by default.}
\item{k}{An integer scalar specifying the number of nearest neighbors to consider when identifying mutual nearest neighbors.}
\item{sigma}{A numeric scalar specifying the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell.}
\item{cos.norm}{A logical scalar indicating whether cosine normalization should be performed.}
\item{svd.dim}{An integer scalar specifying the number of dimensions to use for summarizing biological substructure within each batch.}
\item{order}{An integer vector specifying the order in which batches are to be corrected.}
}

\value{
A named list containing two components:
\describe{
\item{corrected}{A list of length equal to the number of batches, containing matrices of corrected expression values for each cell in each batch.
The order of batches is the same as supplied in \code{...}, and the order of cells in each matrix is also unchanged.
}
\item{mnn}{An integer matrix with two columns, specifying the number of MNNs used for correction of each batch.
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

\section{Further options}{
The input expression values should generally be log-transformed, e.g., log-counts, see \code{\link{normalize}} for details.
By default, a further cosine normalization step is performed on the supplied expression data prior to identifying MNNs, etc.
This can be turned off with \code{cos.norm=FALSE}, though we do not advise doing so.

The function depends on a shared biological manifold, i.e., one or more cell types/states being present in multiple batches.
If this is not true, MNNs may be incorrectly identified.
We protect against this by raising a warning if there are no shared biological subspaces between batches.
We also remove components of the correction vectors that are parallel to the biological subspaces in each batch, to avoid over-correction of interesting biology.
The biological subspace in each batch is identified with SVD, analogous to taking the first several principal components after PCA.
The number of dimensions of this subspace can be controlled with \code{svd.dim}.

Users should note that the order in which batches are corrected will affect the final results.
The first batch in \code{order} is used as the reference batch against which the second batch is corrected.
Corrected values of the second batch are added to the reference batch, against which the third batch is corrected, and so on.
This strategy maximizes the chance of detecting sufficient MNN pairs for stable calculation of correction vectors.
We would consider 20 cells involved in MNN pairs to be the minimum number required for batch correction.
}

\author{
Laleh Haghverdi,
with modifications by Aaron Lun
}

\seealso{
\code{\link[FNN]{get.knnx}}
}

\examples{
B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
out <- mnnCorrect(B1, B2) # corrected values
}
