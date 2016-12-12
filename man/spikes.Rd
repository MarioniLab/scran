\name{Get spikes}
\alias{isSpike<-}
\alias{isSpike<-,SCESet-method}

\title{Set spike-in rows}
\description{Specify the rows in the SCESet corresponding to spike-in transcripts.}

\usage{
\S4method{isSpike}{SCESet}(x) <- value
}

\arguments{
\item{x}{A SCESet object with spike-in data in the \code{colData}.}
\item{value}{
A character vector specifying which control sets are spike-ins.
Alternatively a \code{NULL} value, to remove existing spike-in specifications.
}
}

\details{
This function is deprecated, use \code{\link[scater]{setSpike}} instead.
}

\value{
For \code{isSpike<-}, \code{x} is modified to store a spike-specifying vector in \code{fData(x)$is_feature_spike}.
A logical vector indicating which controls are spike-ins is also stored in the \code{featureControlInfo} slot of \code{x}.
}

\seealso{
\code{\link[scater]{setSpike}},
}

\author{
Aaron Lun
}

\examples{
example(computeSpikeFactors) 
setSpike(y) <- "Spike"
}
