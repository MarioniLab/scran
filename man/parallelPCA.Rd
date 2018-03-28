\name{Parallel analysis}
\alias{parallelPCA}
\alias{parallelPCA,ANY-method}
\alias{parallelPCA,SingleCellExperiment-method}

\title{Parallel analysis for PCA}
\description{Perform a parallel analysis to choose the number of principal components.}

\usage{
\S4method{parallelPCA}{ANY}(x, subset.row=NULL, scale=NULL, value=c("pca", "n", "lowrank"), 
    min.rank=5, max.rank=100, niters=20, keep.perm=FALSE, approximate=FALSE, 
    irlba.args=list(), BPPARAM=SerialParam())

\S4method{parallelPCA}{SingleCellExperiment}(x, ..., subset.row=NULL, 
    value=c("pca", "n", "lowrank"), assay.type="logcounts", 
    get.spikes=FALSE, sce.out=TRUE)
}

\arguments{
\item{x}{A numeric matrix of log-expression values for \code{parallelPCA,ANY-method}, 
or a SingleCellExperiment object containing such values for \code{parallelPCA,SingleCellExperiment-method}.}
\item{subset.row}{A logical, integer or character vector indicating the rows of \code{x} to use for PCA.
All genes are used by default.}
\item{scale}{A numeric vector specifying the scaling to apply to each row of \code{x}, if any.}
\item{value}{A string specifying the type of value to return; the PCs, the number of retained components, or a low-rank approximation.}
\item{min.rank, max.rank}{Integer scalars specifying the minimum and maximum number of PCs to retain.}
\item{niters}{Integer scalar specifying the number of iterations to use for the parallel analysis.}
\item{keep.perm}{Logical scalar indicating whether permutation results should be returned.}
\item{approximate}{A logical scalar indicating whether approximate SVD should be performed via \code{\link[irlba]{irlba}}.}
\item{irlba.args}{A named list of additional arguments to pass to \code{\link[irlba]{irlba}} when \code{approximate=TRUE}.}
\item{BPPARAM}{A BiocParallelParam object.}
\item{...}{Further arguments to pass to \code{denoisePCA,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{get.spikes}{A logical scalar specifying whether spike-in transcripts should be used.
This will be intersected with \code{subset.row} if the latter is specified.}
\item{sce.out}{A logical scalar specifying whether a modified SingleCellExperiment object should be returned.}
}

\value{
For \code{parallelPCA,ANY-method}, a numeric matrix is returned containing the selected PCs (columns) for all cells (rows) if \code{value="pca"}.
If \code{value="n"}, it will return an integer scalar specifying the number of retained components.
If \code{value="lowrank"}, it will return a low-rank approximation of \code{x} with the \emph{same} dimensions.

For \code{parallelPCA,SingleCellExperiment-method}, the return value is the same as \code{parallelPCA,ANY-method} if \code{sce.out=FALSE} or \code{value="n"}.
Otherwise, a SingleCellExperiment object is returned that is a modified version of \code{x}.
If \code{value="pca"}, the modified object will contain the PCs as the \code{"PCA"} entry in the \code{reducedDims} slot.
If \code{value="lowrank"}, it will return a low-rank approximation in \code{assays} slot, named \code{"lowrank"}.

In all cases, the fractions of variance explained by the first \code{max.rank} PCs will be stored as the \code{"percentVar"} attribute in the return value.
If \code{keep.perm=TRUE}, fractions of variance explained by these PCs after each permutation iteration is also recorded as a matrix in \code{"permuted.percentVar"}.
}

\details{
This function performs Horn's parallel analysis to decide how many PCs to retain in a principal components analysis. 
Parallel analysis involves permuting the expression vector for each gene and repeating the PCA to obtain the fractions of variance explained under a random null hypothesis.
The number of PCs to retain is determined by the intersection of the \dQuote{fraction explained} lines on a scree plot.
This is justified as discarding PCs that explain less variance than would be expected under a random model.

This function can be sped up by specifying \code{approximate=TRUE}, which will use approximate strategies for performing the PCA.
Another option is to set \code{BPPARAM} to perform the iterations in parallel.
}

\section{Further comments on parallel analysis}{
Unlike other implementations, this function does \emph{not} compute p-values for the \dQuote{significance} of each principal component.
Such p-values only have a clear interpretation for the first eigenvalue, as the null hypothesis is that of global randomness.
It is not obvious how to perform correction for multiple testing across multiple principal components.
In addition, the number of PCs will depend on the power, which is determined by the number of iterations.

Here, we only use the average of the fractions explained by each successive PC across \code{niters} permutations.
The number of PCs is defined as the last consecutive PC (starting from the first) where the fraction explained is greater than the permuted average.
This is stable to the number of iterations as it uses the expectation of the fractions.
Note that the number of PCs will be coerced to lie between \code{min.rank} and \code{max.rank}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{denoisePCA}}
}

\examples{
# Mocking up some data.
ngenes <- 1000
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, 
            size=1/dispersions), ncol=nsamples)

# Choosing the number of PCs
lcounts <- log2(counts + 1)
parallelPCA(lcounts, min.rank=0, value="n")
}

\references{
Buja A and Eyuboglu N (1992). 
Remarks on Parallel Analysis. 
\emph{Multivariate Behav. Res.}, 27:509-40.
}
