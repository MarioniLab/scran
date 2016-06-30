\name{Spike-in normalization}
\alias{computeSpikeFactors}
\alias{computeSpikeFactors,SCESet-method}

\title{Normalization with spike-in counts}
\description{Compute size factors based on the coverage of spike-in transcripts.}

\usage{
\S4method{computeSpikeFactors}{SCESet}(x, type=NULL, sf.out=FALSE)
}

\arguments{
\item{x}{A SCESet object containing rows corresponding spike-in transcripts.}
\item{type}{A string specifying which spike-in sets to use.}
\item{sf.out}{A logical scalar indicating whether only size factors should be returned.}
}

\details{
The size factor for each cell is defined as the sum of all spike-in counts in each cell.
This is equivalent to normalizing to equalize spike-in coverage between cells.
Size factors are scaled so that the mean of all size factors is unity, for standardization purposes if one were to compare different sets of size factors.

Spike-in counts are assumed to be stored in the rows specified by \code{fData(x)$is_feature_spike} -- see \code{?\link{isSpike<-}} for more details.
If multiple spike-in sets are available, all of them will be used to compute the size factors, if they were all used to set \code{fData(x)$is_feature_spike}.
Individual spike-in sets can be used by specifying the name of the spike-in set in \code{type}.
}

\value{
If \code{sf.out=TRUE}, an object of class \code{x} is returned, containing a numeric vector of size factors for all cells in \code{sizeFactors(x)}.
Otherwise, the size factors are returned directly.
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
