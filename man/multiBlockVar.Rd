\name{multiBlockVar}
\alias{multiBlockVar}

\title{Per-block variance statistics}
\description{Fit a mean-dependent trend to the per-gene variances for each blocking level, and decompose them to biological and technical components.

This has now been deprecated in favor of the \code{block} argument in \code{\link{modelGeneVar}} and related functions.
}

\usage{
multiBlockVar(x, block, make.tech.trend=FALSE, trend.args=list(), 
    dec.args=list(), assay.type="logcounts", ...)
}

\arguments{
\item{x}{A SingleCellExperiment object containing log-normalized expression values, computed with \code{\link{multiBlockNorm}}.}
\item{block}{A factor specifying the blocking level for each cell in \code{x}.}
\item{make.tech.trend}{A logical scalar indicating whether to use \code{\link{makeTechTrend}} to create the mean-variance trend.}
\item{trend.args}{A list of named arguments to pass to \code{\link{trendVar}} if \code{make.tech.trend=FALSE} or \code{\link{makeTechTrend}} otherwise.}
\item{dec.args}{A list of named arguments to pass to \code{\link{decomposeVar}}.}
\item{assay.type}{A string or integer scalar specifying the assay in \code{x} to use for all calculations (except for \code{\link{makeTechTrend}}).}
\item{...}{Additional arguments to pass to \code{\link{combineVar}}.}
}

\details{
This function models the variance of expression in each level of \code{block} separately.
Each subset of cells is passed to \code{\link{trendVar}} (or \code{\link{makeTechTrend}}) to fit a block-specific trend, 
and then passed to \code{\link{decomposeVar}} to obtain block-specific biological and technical components.
Results are consolidated across blocks using the \code{\link{combineVar}} function.
The aim is to enable users to handle differences in the mean-variance relationship across, e.g., different experimental batches.

We assume that the size factors for the endogenous genes have the same mean as the size factors for the spike-ins \emph{within} each block.
This ensures that the spike-in normalized values are comparable to those of the endogenous genes.
Centering should be performed by running \code{\link{multiBlockNorm}} before calling this function.
Otherwise, a warning will be raised about non-centered size factors.
}

\value{
A DataFrame is returned containing all components returned by \code{\link{combineVar}}, in addition to a \code{per.block} column.
This additional column is a DataFrame containing nested DataFrames, each containing a result of \code{\link{decomposeVar}} for the corresponding level of \code{block}.
The trend function from \code{\link{trendVar}} or \code{\link{makeTechTrend}} is also stored as \code{trend} in the metadata of the per-block nested DataFrames.
}

\seealso{
\code{\link{trendVar}},
\code{\link{decomposeVar}},
\code{\link{combineVar}},
\code{\link{multiBlockNorm}}
}

\author{
Aaron Lun
}

\examples{
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.

# Normalizing (gene-based factors for genes, spike-in factors for spike-ins)
y <- computeSumFactors(y) 
y <- computeSpikeFactors(y, general.use=FALSE)

# Setting up the blocking levels.
block <- sample(3, ncol(y), replace=TRUE)
y <- multiBlockNorm(y, block)
out <- multiBlockVar(y, block=block)

# Creating block-level plots.
par(mfrow=c(1,3))
is.spike <- isSpike(y)
for (x in as.character(1:3)) {
    current <- out$per.block[[x]]
    plot(current$mean, current$total, col="black", pch=16)
    points(current$mean[is.spike], current$total[is.spike], col="red", pch=16)
    curve(metadata(current)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
}
}

\references{
Lun ATL, McCarthy DJ and Marioni JC (2016). 
A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.
\emph{F1000Res.} 5:2122
}

\keyword{variance}

