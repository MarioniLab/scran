\name{Parallel analysis}
\alias{parallelPCA}
\alias{parallelPCA,ANY-method}
\alias{parallelPCA,SingleCellExperiment-method}

\title{Parallel analysis for PCA}
\description{Perform a parallel analysis to choose the number of principal components.}

\usage{
\S4method{parallelPCA}{ANY}(x, subset.row=NULL, value=c("pca", "n", "lowrank"), 
    min.rank=5, max.rank=100, niters=50, threshold=0.1, BSPARAM=ExactParam(), 
    BPPARAM=SerialParam())

\S4method{parallelPCA}{SingleCellExperiment}(x, ..., subset.row=NULL, 
    value=c("pca", "n", "lowrank"), assay.type="logcounts", 
    get.spikes=FALSE, sce.out=TRUE)
}

\arguments{
\item{x}{A numeric matrix of log-expression values for \code{parallelPCA,ANY-method}, 
or a SingleCellExperiment object containing such values for \code{parallelPCA,SingleCellExperiment-method}.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{value}{A string specifying the type of value to return; the PCs, the number of retained components, or a low-rank approximation.}
\item{min.rank, max.rank}{Integer scalars specifying the minimum and maximum number of PCs to retain.}
\item{niters}{Integer scalar specifying the number of iterations to use for the parallel analysis.}
\item{threshold}{Numeric scalar representing the \dQuote{p-value} threshold above which PCs are to be ignored.}
\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how the iterations should be paralellized.}
\item{...}{Further arguments to pass to \code{denoisePCA,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{get.spikes}{See \code{?"\link{scran-gene-selection}"}.}
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
Fractions of variance explained by these PCs after each permutation iteration are also recorded as a matrix in \code{"permuted.percentVar"}.
}

\details{
This function performs Horn's parallel analysis to decide how many PCs to retain in a principal components analysis. 
Parallel analysis involves permuting the expression vector for each gene and repeating the PCA to obtain the fractions of variance explained under a random null model. 
The number of PCs to retain is determined by the intersection of the \dQuote{fraction explained} lines on a scree plot.
This is justified as discarding PCs that explain less variance than would be expected under a random model.

In practice, we discard all PCs from the first PC that has a fraction explained \emph{similar} to that under the null.
A PC is considered similar if the permuted fractions exceed the observed fraction in more than \code{threshold} of iterations.
(For want of a better word, we have described this as a \dQuote{p-value} threshold, though it is not interpretable as a measure of significance.)
This is a more conservative criterion than discarding PCs with fractions below the average null fraction, which tends to overstate the rank in noisy datasets.
Note that the number of PCs will be coerced to lie between \code{min.rank} and \code{max.rank}.

This function can be sped up by specifying \code{BSPARAM=IrlbaParam()} or similar, 
to use approximate strategies for performing the PCA.
Another option is to set \code{BPPARAM} to perform the iterations in parallel.
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
