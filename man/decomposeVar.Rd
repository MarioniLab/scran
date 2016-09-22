\name{decomposeVar}
\alias{decomposeVar}
\alias{decomposeVar,matrix,list-method}
\alias{decomposeVar,SCESet,list-method}

\title{Decompose the gene-level variance}
\description{Decompose the gene-specific variance into biological and technical components for single-cell RNA-seq data.}

\usage{
\S4method{decomposeVar}{matrix,list}(x, fit, design=NA, subset.row=NULL)
\S4method{decomposeVar}{SCESet,list}(x, fit, subset.row=NULL, ..., assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{x}{
    A numeric matrix of normalized log-expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SCESet object containing such a matrix.
}
\item{fit}{
    A list containing the output of \code{\link{trendVar}}, run on log-expression values for spike-in genes.
}
\item{design}{A numeric matrix describing the systematic factors contributing to expression in each cell.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to \code{decomposeVar,matrix,list-method}.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\details{
This function computes the variance of the log-CPMs for each endogenous gene.
The technical component of the variance for each gene is determined by interpolating the fitted trend in \code{fit} at the mean log-CPM for that gene.
This represents variance due to sequencing noise, variability in capture efficiency, etc.
The biological component is determined by subtracting the technical component from the total variance.

Highly variable genes (HVGs) can be identified as those with large biological components.
Unlike other methods for decomposition, this approach estimates the variance of the log-CPMs rather than of the counts themselves.
The log-transformation blunts the impact of large positive outliers and ensures that the HVG list is not dominated by outliers.
Interpretation is not compromised -- HVGs will still be so, regardless of whether counts or log-CPMs are considered.

The design matrix can be set if there are factors that should be blocked, e.g., batch effects, known (and uninteresting) clusters.
If \code{NULL}, it will be set to an all-ones matrix, i.e., all cells are replicates.
If \code{NA}, it will be extracted from \code{fit$design}, assuming that the same cells were used to fit the trend.

Users can also directly specify which rows to use with \code{subset.row}.
This is equivalent to running \code{decomposeVar} on \code{x[subset.row,]}, but is more efficient as it avoids the construction of large temporary matrices.
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
set.seed(100)

nspikes <- ncells <- 200
spike.means <- 2^runif(nspikes, 3, 8)
spike.disp <- 100/spike.means + 0.5
spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)

ngenes <- 10000
cell.means <- 2^runif(ngenes, 2, 10)
cell.disp <- 100/cell.means + 0.5
cell.data <- matrix(rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)

combined <- rbind(cell.data, spike.data)
colnames(combined) <- seq_len(ncells)
rownames(combined) <- seq_len(nrow(combined))
y <- newSCESet(countData=combined)
y <- calculateQCMetrics(y, list(Spike=rep(c(FALSE, TRUE), c(ngenes, nspikes))))
isSpike(y) <- "Spike"

# Normalizing.
y <- computeSumFactors(y)
y <- computeSpikeFactors(y, general.use=FALSE)
y <- normalize(y)

# Decomposing technical and biological noise.
fit <- trendVar(y)
results <- decomposeVar(y, fit)

plot(results$mean, results$total)
o <- order(results$mean)
lines(results$mean[o], results$tech[o], col="red", lwd=2)

plot(results$mean, results$bio)
}

\keyword{variance}

