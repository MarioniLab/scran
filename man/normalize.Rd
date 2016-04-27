\name{Normalized Expression}
\alias{normalize}
\alias{normalize,SCESet-method}

\title{Compute normalized expression values}
\description{Compute (log-)expression values from read counts, using pre-calculated size factors for normalization.}

\usage{
\S4method{normalize}{SCESet}(object, log=TRUE, prior.count=1, separate.spikes=TRUE)
}

\arguments{
\item{object}{A matrix of read counts, or a SCESet object with an assay named \code{"counts"}.}
\item{log}{A logical scalar specifying whether the expression should be log-transformed.}
\item{prior.count}{A numeric scalar indicating the prior count to add prior to log-transformation, to avoid undefined values from zero counts.}
\item{separate.spikes}{A logical scalar indicating whether spike-in counts should be normalized separately.}
}

\details{
This function computes normalized log-expression values by adding \code{prior.count} to each count, dividing by the \code{size.factor} for that cell, and log-transforming.
Size factors are taken from the appropriate field in the \code{colData} of \code{object}.
These size factors can be computed with a number of functions like \code{\link{computeSumFactors}} or \code{\link{computeSpikeFactors}}.

If spike-in counts are present in the SCESet object, these will also be converted into normalized values.
If \code{separate.spikes=FALSE}, this is done with the same set of size factors that was used for the endogenous genes.
Otherwise, a separate set of spike-in size factors will be used instead -- these are defined by calling \code{\link{computeSpikeFactors}}.

All size factors are mean-centered so that their geometric mean is equal to unity prior to computing normalized expression values.
This ensures that expression values are roughly comparable when different sets of size factors are used.
}

\section{Why spikes can be normalized separately}{
In most cases, it does not make sense to normalize spike-in counts with size factors computed from endogenous genes.
This is because the spike-in counts do not (generally) depend on the total amount of endogenous RNA, whereas the size factors do.
As such, normalizing the former with the latter would be inappropriate -- cells with a lot of endogenous RNA would scale down the spike-in counts, even if the same amount of spike-in RNA was added, captured and sequenced in each cell.

Instead, normalization of the spike-in counts should be performed using size factors computed from those counts, i.e., with \code{\link{computeSpikeFactors}}.
This is the default setting when \code{separate.spikes=TRUE}.
Normalized log-expression values are made to be roughly comparable to those of endogenous genes, by ensuring \emph{all} sets of size factors are mean-centered prior to normalization.
}

\value{
For \code{normalize,SCESet-method}, a SCESet object is returned with an additional assay named \code{"exprs"}.
This contains normalized log-expression values for the endogenous genes.
If spike-ins are present, normalized values for the spike-in transcripts are stored in the \code{norm.spikes} field of the \code{colData}.
}

\seealso{
\code{\link{cpm}},
\code{\link{SCESet}},
\code{\link{computeSumFactors}},
\code{\link{computeSpikeFactors}}
}

\author{
Aaron Lun
}

\examples{
set.seed(100)
popsize <- 10
ngenes <- 1000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=10*all.facs, size=1), ncol=popsize, byrow=TRUE)
spikes <- matrix(rnbinom(100*popsize, mu=10*all.facs, size=0.5), ncol=popsize, byrow=TRUE)

combined <- rbind(counts, spikes)
colnames(combined) <- seq_len(popsize)
rownames(combined) <- seq_len(nrow(combined))
y <- newSCESet(countData=combined)
isSpike(y) <- rep(c(FALSE, TRUE), c(ngenes, 100))

sizeFactors(y) <- colSums(combined) # Library size normalization, basically.
y <- normalize(y)
exprs(y)[1:10,]

y <- computeSpikeFactors(y)
y <- normalize(y)
exprs(y)[1:10,]
}
