\name{decomposeVar}
\alias{decomposeVar}
\alias{decomposeVar,matrix,list-method}
\alias{decomposeVar,SCESet,list-method}

\title{Decompose the gene-level variance}
\description{Decompose the gene-specific variance into biological and technical components for single-cell RNA-seq data.}

\usage{
\S4method{decomposeVar}{matrix,list}(x, fit, design=NA)
\S4method{decomposeVar}{SCESet,list}(x, fit, ..., assay="exprs", get.spikes=FALSE)
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
}

\value{
A data frame is returned where each row corresponds to and is named after a row of \code{x}.
This contains the numeric fields:
\describe{
\item{\code{mean}:}{Mean normalized log-count per gene.}
\item{\code{total}:}{Variance of the normalized log-counts per gene.}
\item{\code{bio}:}{Biological component of the variance.}
\item{\code{tech}:}{Technical component of the variance.}
}
Rows corresponding to spike-in transcripts are set to \code{NA} unless \code{get.spikes=TRUE}. 
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

nspikes <- ncells <- 100
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
isSpike(y) <- rep(c(FALSE, TRUE), c(ngenes, nspikes))

# Normalizing.
y <- computeSpikeFactors(y) # or computeSumFactors
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

% # Compare the pair - 'a' has higher CV2 and lower var(log), while 'b' has lower CV2 and higher var(log).
% a <- integer(100)
% a[1] <- 10000
% mean(a)
% var(a)/mean(a)^2
% var(log(a+1))
% b <- rep(c(0, 200), each=50)
% mean(b)
% var(b)/mean(b)^2
% var(log(b+1))
%
% # Another way of looking at it, assuming an outlier of size 'x' and n<-100 cells.
% # Variance of the log for NB-distributed (disp=0.5) spike-ins is greater than the outlier at each mean.
% # The same cannot be said for the CV2.
% n <- 100
% par(mfrow=c(1,2))
% curve(log(1+1/(x/n)+0.5), xlim=c(1, 1000), ylim=c(0,2)) # variance of the log (based on second moment of log-normal -> NB)
% curve(log(x)^2/100, add=TRUE, col="red") # variance of log for one outlier 'x'; analytically derived from mean 'x/n' + Bessel's correction.
% curve(1/(x/n) + 0.5, xlim=c(1, 1000)) # CV2 of the spike-in counts.
% abline(h=n, col="blue") # CV2; divide variance by '(x/n)^2'.

