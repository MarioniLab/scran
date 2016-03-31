\name{Spike-in normalization}
\alias{computeSpikeFactors}
\alias{computeSpikeFactors,SCESet-method}

\title{Normalization with spike-in counts}
\description{Compute size factors based on the coverage of spike-in transcripts.}

\usage{
\S4method{computeSpikeFactors}{SCESet}(x)
}

\arguments{
\item{x}{A SCESet object containing rows corresponding spike-in transcripts.}
}

\details{
The size factor for each cell is defined as the sum of all spike-in counts in each cell.
This is equivalent to normalizing to equalize spike-in coverage between cells.
Spike-in counts are assumed to be stored in rows with \code{fData(x)$is_control_feature} -- see \code{?\link{isSpike<-}} for more details.
Note that the geometric mean of all size factors is set to unity, for standardization purposes if one were to compare different sets of size factors.
}

\value{
An object of class \code{x} is returned, containing a numeric vector of size factors for all cells in \code{sizeFactors(x)}.
}

\author{
Aaron Lun
}

\examples{
set.seed(100)
popsize <- 200
ngenes <- 1000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs*10, size=1), ncol=popsize, byrow=TRUE)
spikes <- matrix(rnbinom(100*popsize, mu=all.facs*10, size=0.5), ncol=popsize, byrow=TRUE)

combined <- rbind(counts, spikes)
colnames(combined) <- seq_len(popsize)
rownames(combined) <- seq_len(nrow(combined))
y <- newSCESet(countData=combined)
isSpike(y) <- rep(c(FALSE, TRUE), c(ngenes, 100))
out.facs <- computeSpikeFactors(y)
}

\seealso{
\code{\link{SCESet}}
}

\keyword{normalization}
