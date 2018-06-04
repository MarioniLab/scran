\name{multiBatchNorm}
\alias{multiBatchNorm}

\title{Per-batch scaling normalization}
\description{Perform scaling normalization within each batch to provide comparable results to the lowest-coverage batch.}

\usage{
multiBatchNorm(..., assay.type="counts", norm.args=list(), 
    subset.row=NULL, min.mean=1)
}

\arguments{
\item{...}{Two or more SingleCellExperiment objects containing counts and size factors.
Each object is assumed to represent one batch.}
\item{assay.type}{A string specifying which assay values contains the counts.}
\item{norm.args}{A named list of further arguments to pass to \code{\link[scater]{normalize}}.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{min.mean}{A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.}
}

\details{
When performing integrative analyses of multiple batches, it is often the case that different batches have large differences in coverage.
This function removes systematic differences in coverage across batches to simplify downstream comparisons.
It does so by resaling the size factors using median-based normalization on the ratio of the average counts between batches.
This is equivalent to the between-cluster normalization performed in \code{\link{computeSumFactors}}.

This function will adjust the size factors so that counts in high-coverage batches are scaled \emph{downwards} to match the coverage of the most shallow batch.
The \code{\link[scater]{normalize}} function will then add the same pseudo-count to all batches before log-transformation.
By scaling downwards, we favour stronger squeezing of log-fold changes from the pseudo-count, mitigating any differences in variance between batches.
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

\keyword{normalization}
