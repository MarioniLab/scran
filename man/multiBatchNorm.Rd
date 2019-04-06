\name{multiBatchNorm}
\alias{multiBatchNorm}

\title{Per-batch scaling normalization}
\description{Perform scaling normalization within each batch to provide comparable results to the lowest-coverage batch.}

\usage{
multiBatchNorm(..., assay.type="counts", norm.args=list(), min.mean=1,
    subset.row=NULL)
}

\arguments{
\item{...}{Two or more SingleCellExperiment objects containing counts and size factors.
Each object is assumed to represent one batch.}
\item{assay.type}{A string specifying which assay values contains the counts.}
\item{norm.args}{A named list of further arguments to pass to \code{\link[scater]{normalize}}.}
\item{min.mean}{A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
}

\details{
When performing integrative analyses of multiple batches, it is often the case that different batches have large differences in coverage.
This function removes systematic differences in coverage across batches to simplify downstream comparisons.
It does so by resaling the size factors using median-based normalization on the ratio of the average counts between batches.
This is roughly equivalent to the between-cluster normalization performed in \code{\link{computeSumFactors}}.

This function will adjust the size factors so that counts in high-coverage batches are scaled \emph{downwards} to match the coverage of the most shallow batch.
The \code{\link[scater]{normalize}} function will then add the same pseudo-count to all batches before log-transformation.
By scaling downwards, we favour stronger squeezing of log-fold changes from the pseudo-count, mitigating any technical differences in variance between batches.
Of course, genuine biological differences will also be shrunk, but this is less of an issue for upregulated genes with large counts.

Running this function is preferred over running \code{\link{normalize}} directly when computing log-normalized values for use in \code{\link{mnnCorrect}} or \code{\link{fastMNN}}.
This is because, in most cases, size factors will be computed within each batch;
their direct application in \code{\link{normalize}} will not account for scaling differences between batches.
In contrast, \code{multiBatchNorm} will rescale the size factors so that they are comparable across batches.

If spike-in transcripts are present, these should be the same across all batches.
Spike-in size factors are rescaled separately from those of the endogenous genes, to reflect differences in spike-in quantities across batches.
Conversely, spike-in transcripts are not used to compute the rescaling factors for endogenous genes.

Only genes with library size-adjusted average counts greater than \code{min.mean} will be used for computing the rescaling factors.
This improves precision and avoids problems with discreteness.
Users can also set \code{subset.row} to restrict the set of genes used for computing the rescaling factors.
However, this only affects the rescaling of the size factors - normalized values for \emph{all} genes will still be returned.
}

\value{
A list of SingleCellExperiment objects with normalized log-expression values in the \code{"logcounts"} assay (depending on values in \code{norm.args}).
}

\author{
Aaron Lun
}

\seealso{
\code{\link[scater]{normalize}},
\code{\link{mnnCorrect}},
\code{\link{fastMNN}}
}

\examples{
\dontrun{
d1 <- matrix(rnbinom(50000, mu=10, size=1), ncol=100)
sce1 <- SingleCellExperiment(list(counts=d1))
sizeFactors(sce1) <- runif(ncol(d1))

d2 <- matrix(rnbinom(20000, mu=50, size=1), ncol=40)
sce2 <- SingleCellExperiment(list(counts=d2))
sizeFactors(sce2) <- runif(ncol(d2))

out <- multiBatchNorm(sce1, sce2)
summary(sizeFactors(out[[1]]))
summary(sizeFactors(out[[2]]))
}
}

\keyword{normalization}
