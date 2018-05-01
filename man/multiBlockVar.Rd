\name{multiBlockVar}
\alias{multiBlockVar}

\title{Per-block variance statistics}
\description{Fit a mean-dependent trend to the per-gene variances for each blocking level, and decompose them to biological and technical components.}

\usage{
multiBlockVar(x, block, trend.args=list(), dec.args=list(), ...)
}

\arguments{
\item{x}{A SingleCellExperiment object containing log-normalized expression values,
computed from size factors centred in each block - see \code{\link[scater]{normalize}} for details.
}
\item{block}{A factor specifying the blocking level for each cell.}
\item{trend.args}{A list of named arguments to pass to \code{\link{trendVar}}.}
\item{dec.args}{A list of named arguments to pass to \code{\link{decomposeVar}}.}
\item{...}{Additional arguments to pass to \code{\link{combineVar}}.}
}

\details{
This function models the variance of expression in each level of \code{block} separately.
Each subset of cells is passed to \code{\link{trendVar}} to fit a block-specific trend, 
and then passed to \code{\link{decomposeVar}} to obtain block-specific biological and technical components.
Results are consolidated across blocks using the \code{\link{combineVar}} function.
The aim is to enable users to handle differences in the mean-variance relationship across, e.g., different experimental batches.

We assume that the size factors are centred \emph{within} each block when computing log-normalized expression values.
This preserves the scale of the counts within each block, and ensures that the spike-in normalized values are comparable to those of the endogenous genes.
Centring can be performed by setting the \code{size_factor_grouping} argument in \code{\link[scater]{normalize}}.
Otherwise, a warning will be raised about non-centred size factors.
}

\value{
A DataFrame is returned containing all components returned by \code{\link{combineVar}}, in addition to a \code{per.block} column.
This additional column is a DataFrame containing nested DataFrames, each containing a result of \code{\link{decomposeVar}} for the corresponding level of \code{block}.
The trend function from \code{\link{trendVar}} is also stored as \code{trend} in the metadata of the per-block nested DataFrames.
}

\seealso{
\code{\link{trendVar}},
\code{\link{decomposeVar}},
\code{\link{combineVar}},
\code{\link[scater]{normalize}}
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
y <- normalize(y, size_factor_grouping=block)
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

