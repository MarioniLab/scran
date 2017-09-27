\name{Denoise with PCA}
\alias{denoisePCA}
\alias{denoisePCA,ANY-method}
\alias{denoisePCA,SingleCellExperiment-method}

\title{Denoise expression with PCA}
\description{Denoise log-expression data by removing principal components corresponding to technical noise.}

\usage{
\S4method{denoisePCA}{ANY}(x, technical, design=NULL, subset.row=NULL,
    value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100, 
    preserve.dim=FALSE, approximate=FALSE, rand.seed=1000)

\S4method{denoisePCA}{SingleCellExperiment}(x, ..., subset.row=NULL, 
    value=c("pca", "n", "lowrank"), assay.type="logcounts", 
    get.spikes=FALSE, sce.out=TRUE)
}

\arguments{
\item{x}{A numeric matrix of log-expression values for \code{denoisePCA,ANY-method}, 
or a SingleCellExperiment object containing such values for \code{denoisePCA,SingleCellExperiment-method}.}
\item{technical}{A function that computes the technical component of the variance for a gene with a given mean (log-)expression, see \code{?\link{trendVar}}.}
\item{design}{A numeric matrix containing the experimental design.
If \code{NULL}, all cells are assumed to belong to a single group.}
\item{subset.row}{A logical, integer or character vector indicating the rows of \code{x} to use.
All genes are used by default.}
\item{value}{A string specifying the type of value to return; the PCs, the number of retained components, or a low-rank approximation.}
\item{min.rank, max.rank}{Integer scalars specifying the minimum and maximum number of PCs to retain.}
\item{preserve.dim}{A logical scalar indicating whether the dimensions should be preserved when \code{subset.row} is not \code{NULL}.}
\item{approximate}{A logical scalar indicating whether approximate SVD should be performed via \code{\link[irlba]{irlba}}.}
\item{rand.seed}{A numeric scalar specifying the seed for approximate PCA when \code{approximate=TRUE}.
This can be set to \code{NA} to use the existing session seed.}
\item{...}{Further arguments to pass to \code{denoisePCA,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{get.spikes}{A logical scalar specifying whether spike-in transcripts should be used.
This will be intersected with \code{subset.row} if the latter is specified.}
\item{sce.out}{A logical scalar specifying whether a modified SingleCellExperiment object should be returned.}
}

\value{
For \code{denoisePCA,ANY-method}, a numeric matrix is returned containing the selected PCs (columns) for all cells (rows) if \code{value="pca"}.
If \code{value="n"}, it will return an integer scalar specifying the number of retained components.
If \code{value="lowrank"}, it will return a low-rank approximation of \code{x} with the same dimensions.

For \code{denoisePCA,SingleCellExperiment-method}, the return value is the same as \code{denoisePCA,ANY-method} if \code{sce.out=TRUE} or \code{value="n"}.
Otherwise, a SingleCellExperiment object is returned that is a modified version of \code{x}.
If \code{value="pca"}, the modified object will contain the PCs as the \code{"PCA"} entry in the \code{reducedDims} slot.
If \code{value="lowrank"}, it will return a low-rank approximation in \code{assays} slot, named \code{"lowrank"}.
}

\details{
This function performs a principal components analysis to reduce random technical noise in the data.
Random noise is uncorrelated across genes and should be captured by later PCs, as the variance in the data explained by any single gene is low.
In contrast, biological substructure should be correlated and captured by earlier PCs, as this explains more variance for sets of genes.
The idea is to discard later PCs to remove technical noise and improve the resolution of substructure.

The choice of the number of PCs to discard is based on the estimates of technical variance in \code{technical}.
This uses the trend function obtained from \code{\link{trendVar}} to compute the technical component for each gene, based on its mean abundance.
The overall technical variance is estimated by summing the values across genes.
Genes with negative biological components are ignored during downstream analyses to ensure that the total variance is greater than the technical variance.

The function works by assuming that the first X PCs contain all of the biological signal, while the remainder contains technical noise.
For a given value of X, a prediction for the total technical variance is calculated, based on the variance explained by the later PCs.
(It also accounts for a non-zero technical component in the variance explained by first X PCs.)
A value of X is found such that the predicted technical variance equals the estimated technical variance.
Note that X will be coerced to lie between \code{min.rank} and \code{max.rank}.

Only the first X PCs are reported if \code{value="pca"}.
If \code{value="n"}, the value of X is directly reported, which avoids computing the PCs if only the rank is desired.
If \code{value="lowrank"}, a low-rank approximation of the original matrix is computed using only the first X components.
This is useful for denoising prior to downstream applications that expect gene-wise expression profiles.

If \code{value="lowrank"}, genes with negative components are still reported but are assigned expression values of zero for all cells.
If \code{subset.row} is not \code{NULL}, genes that are not in the selected set are removed by default.
However, if \code{preserve.dim=TRUE}, the unselected genes will still be reported and are assigned all-zero expression profiles.

If \code{design} is specified, the residuals of a linear model fitted to each gene are computed.
Because variances computed from residuals are usually underestimated, the residuals are scaled up so that their variance is equal to the residual variance of the model fit.
This ensures that the sum of variances is not understated, which would lead to more PCs being discarded than appropriate.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{trendVar}},
\code{\link{decomposeVar}}
}

\examples{
# Mocking up some data.
ngenes <- 1000
is.spike <- 1:100
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), ncol=nsamples)
rownames(counts) <- paste0("Gene", seq_len(ngenes))

# Fitting a trend.
lcounts <- log2(counts + 1)
fit <- trendVar(lcounts, subset.row=is.spike)

# Denoising (not including the spike-ins in the PCA;
# spike-ins are automatically removed with the SingleCellExperiment method). 
pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=-is.spike)
dim(pcs)

# With a design matrix.
design <- model.matrix(~factor(rep(0:1, length.out=nsamples)))
fit3 <- trendVar(lcounts, design=design, subset.row=is.spike)
pcs3 <- denoisePCA(lcounts, technical=fit3$trend, design=design, subset.row=-is.spike)
dim(pcs3)
}

