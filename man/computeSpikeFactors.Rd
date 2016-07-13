\name{Spike-in normalization}
\alias{computeSpikeFactors}
\alias{computeSpikeFactors,SCESet-method}

\title{Normalization with spike-in counts}
\description{Compute size factors based on the coverage of spike-in transcripts.}

\usage{
\S4method{computeSpikeFactors}{SCESet}(x, type=NULL, sf.out=FALSE, general.use=TRUE)
}

\arguments{
\item{x}{A SCESet object containing rows corresponding spike-in transcripts.}
\item{type}{A character vector specifying which spike-in sets to use.}
\item{sf.out}{A logical scalar indicating whether only size factors should be returned.}
\item{general.use}{A logical scalar indicating whether the size factors should be stored for general use by all genes.}
}

\details{
The size factor for each cell is defined as the sum of all spike-in counts in each cell.
This is equivalent to normalizing to equalize spike-in coverage between cells.
Size factors are scaled so that the mean of all size factors is unity, for standardization purposes if one were to compare different sets of size factors.

Spike-in counts are assumed to be stored in the rows specified by \code{isSpike(x)}.
This specification should have been performed by supplying the names of the spike-in sets -- see \code{?\link{isSpike<-}} for more details.
By default, if multiple spike-in sets are available, all of them will be used to compute the size factors.
The function can be restricted to a subset of the spike-ins by specifying the names of the desired spike-in sets in \code{type}.

By default, the function will store several copies of the same size factors in the output object.
One copy will be stored in \code{sizeFactors(x)} for normalization of all genes -- this can be disabled by setting \code{general.use=FALSE}.
One copy will also be stored in \code{sizeFactors(x, type=s)}, where \code{s} is the name of a spike-in set in \code{type}.
(If \code{type=NULL}, a copy is stored for every spike-in set, as all of them would be used to compute the size factors.)
Separate storage allows spike-in-specific normalization in \code{\link[scater]{normalize,SCESet-method}}.
}

\value{
If \code{sf.out=TRUE}, a numeric vector of size factors is returned directly.

Otherwise, an object of class \code{x} is returned, containing size factors for all cells.
A copy of the vector is stored for each spike-in set that was used to compute the size factors.
If \code{general.use=TRUE}, a copy is also stored for use by non-spike-in genes.
}

\author{
Aaron Lun
}

\examples{
# Setting up an example.
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
y <- calculateQCMetrics(y, list(IAmASpike=rep(c(FALSE, TRUE), c(ngenes, 100))))
isSpike(y) <- "IAmASpike"

# 
y <- computeSpikeFactors(y)
sizeFactors(y)
sizeFactors(y, type="IAmASpike")

# general.use=FALSE does not modify general size factors
y2 <- y
sizeFactors(y2) <- 1
sizeFactors(y2, type="IAmASpike") <- 1
y2 <- computeSpikeFactors(y2, general.use=FALSE)
sizeFactors(y2)
sizeFactors(y2, type="IAmASpike")

}

\seealso{
\code{\link{SCESet}}
}

\keyword{normalization}
