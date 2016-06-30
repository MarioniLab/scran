\name{Get spikes}
\alias{spikes}
\alias{spikes,SCESet-method}
\alias{isSpike}
\alias{isSpike,SCESet-method}
\alias{isSpike<-}
\alias{isSpike<-,SCESet,logical-method}
\alias{isSpike<-,SCESet,character-method}
\alias{isSpike<-,SCESet,NULL-method}

\title{Construct the spike-in matrix}
\description{Identify rows in the SCESet corresponding to spike-in transcripts, and retrieve a matrix of counts or normalized expression values for those rows.}

\usage{
\S4method{isSpike}{SCESet}(x, type=NULL)
\S4method{spikes}{SCESet}(x, assay="counts", type=NULL)

\S4method{isSpike}{SCESet,logical}(x) <- value
\S4method{isSpike}{SCESet,character}(x) <- value
\S4method{isSpike}{SCESet,NULL}(x) <- value
}

\arguments{
\item{x}{A SCESet object with spike-in data in the \code{colData}.}
\item{type}{A string specifying which spike-in set should be extracted.}
\item{assay}{A string specifying whether counts or normalized expression values are to be extracted.}
\item{value}{
A logical vector specifying which rows correspond to spike-ins,
or a character vector specifying which control sets are spike-ins.
Alternatively a \code{NULL} value, to remove existing spike-in specifications.
}
}

\details{
The \code{isSpikes} and \code{spikes} methods indicate which rows correspond to spike-ins and their expression values, respectively.
If multiple spike-in sets are available, users can extract information for a specific set by supplying the name of the set in \code{type}.
(By default, rows from all spike-in sets are extracted when \code{type=NULL}.)
If \code{assay="exprs"}, users should have run \code{x} through \code{\link{normalize}}.

Users can set a logical vector in \code{isSpike(x)} indicating whether each row of \code{x} is a spike-in control.
Alternatively, if \code{\link[scater]{calculateQCMetrics}} has been used to process the data, the name(s) of the spike-in set(s) can be supplied.
In the latter case, the logical vector will be automatically constructed and include rows from all supplied sets.
Setting \code{NULL} will clear all existing spike-in information in \code{x}.
}

\value{
For \code{spikes}, a numeric matrix of counts or normalized expression values, with one column per cell and one row per spike-in transcript.

For \code{isSpike}, a logical vector indicating which rows are spike-ins (or \code{NULL}, if this information was not stored in \code{x}).
For \code{isSpike<-}, \code{x} is modified to store a spike-specifying vector in \code{fData(x)$is_feature_spike}.
}

\seealso{
\code{\link{normalize}},
\code{\link[scater]{calculateQCMetrics}},
\code{\link{SCESet}}
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

y <- computeSpikeFactors(y)
y <- normalize(y)
spikes(y)[1:10,]
spikes(y, assay="exprs")[1:10,]
isSpike(y)
}
