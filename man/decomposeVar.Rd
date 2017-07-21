\name{decomposeVar}
\alias{decomposeVar}
\alias{decomposeVar,ANY,list-method}
\alias{decomposeVar,SCESet,list-method}

\title{Decompose the gene-level variance}
\description{Decompose the gene-specific variance into biological and technical components for single-cell RNA-seq data.}

\usage{
\S4method{decomposeVar}{ANY,list}(x, fit, design=NA, subset.row=NULL, ...)

\S4method{decomposeVar}{SCESet,list}(x, fit, subset.row=NULL, ..., assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{x}{
    A numeric matrix-like object of normalized log-expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SCESet object containing such a matrix.
}
\item{fit}{
    A list containing the output of \code{\link{trendVar}}, run on log-expression values for spike-in genes.
}
\item{design}{A numeric matrix describing the uninteresting factors contributing to expression in each cell.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{
For \code{decomposeVar,matrix,list-method}, additional arguments to pass to \code{\link{testVar}}.
For \code{decomposeVar,SCESet,list-method}, additional arguments to pass to the matrix method.
}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\details{
This function computes the variance of the normalized log-counts for each endogenous gene.
The technical component of the variance for each gene is determined by interpolating the fitted trend in \code{fit} at the mean log-count for that gene.
This represents variance due to sequencing noise, variability in capture efficiency, etc.
The biological component is determined by subtracting the technical component from the total variance.

Highly variable genes (HVGs) can be identified as those with large biological components.
Unlike other methods for decomposition, this approach estimates the variance of the log-counts rather than of the counts themselves.
The log-transformation blunts the impact of large positive outliers and ensures that the HVG list is not dominated by outliers.
Interpretation is not compromised -- HVGs will still be so, regardless of whether counts or log-counts are considered.

The design matrix can be set if there are factors that should be blocked, e.g., batch effects, known (and uninteresting) clusters.
If \code{NULL}, it will be set to an all-ones matrix, i.e., all cells are replicates.
If \code{NA}, it will be extracted from \code{fit$design}, assuming that the same cells were used to fit the trend.

Users can also directly specify which rows to use with \code{subset.row}.
This is equivalent to running \code{decomposeVar} on \code{x[subset.row,]}, but is more efficient as it avoids the construction of large temporary matrices.

If \code{assay="exprs"} and the size factors are not centred at unity, a warning will be raised - see \code{?\link{trendVar}} for details.
}

\value{
A data frame is returned where each row corresponds to and is named after a row of \code{x} 
(if \code{subset.row=NULL}; otherwise, each row corresponds to an element of \code{subset.row}).
This contains the numeric fields:
\describe{
\item{\code{mean}:}{Mean normalized log-count per gene.}
\item{\code{total}:}{Variance of the normalized log-counts per gene.}
\item{\code{bio}:}{Biological component of the variance.}
\item{\code{tech}:}{Technical component of the variance.}
\item{\code{p.value, FDR}:}{Raw and adjusted p-values for the test against the null hypothesis that \code{bio=0}.}
}
Rows corresponding to spike-in transcripts have their p-value and FDR fields set to \code{NA} unless \code{get.spikes=TRUE}. 
}

\seealso{
\code{\link{trendVar}},
\code{\link{testVar}}
}

\author{
Aaron Lun
}

\examples{
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.
y <- computeSumFactors(y) # Size factors for the the endogenous genes.
y <- computeSpikeFactors(y, general.use=FALSE) # Size factors for spike-ins. 
y <- normalize(y) # Normalizing the counts by the size factors.

# Decomposing technical and biological noise.
fit <- trendVar(y)
results <- decomposeVar(y, fit)
head(results)

plot(results$mean, results$total)
o <- order(results$mean)
lines(results$mean[o], results$tech[o], col="red", lwd=2)

plot(results$mean, results$bio)

# A trend fitted to endogenous genes can also be used, pending assumptions.
fit.g <- trendVar(y, use.spikes=FALSE)
results.g <- decomposeVar(y, fit.g)
head(results.g)
}

\references{
Lun ATL, McCarthy DJ and Marioni JC (2016). 
A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.
\emph{F1000Res.} 5:2122
}

\keyword{variance}

