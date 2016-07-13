\name{Get spikes}
\alias{spikes}
\alias{spikes,SCESet-method}
\alias{isSpike}
\alias{isSpike,SCESet-method}
\alias{isSpike<-}
\alias{isSpike<-,SCESet,character-method}
\alias{isSpike<-,SCESet,NULL-method}

\title{Construct the spike-in matrix}
\description{Identify rows in the SCESet corresponding to spike-in transcripts, and retrieve a matrix of counts or normalized expression values for those rows.}

\usage{
\S4method{isSpike}{SCESet}(x, type=NULL)
\S4method{spikes}{SCESet}(x, assay="counts", type=NULL)

\S4method{isSpike}{SCESet,character}(x) <- value
\S4method{isSpike}{SCESet,NULL}(x) <- value
}

\arguments{
\item{x}{A SCESet object with spike-in data in the \code{colData}.}
\item{type}{A character vector specifying which spike-in set(s) should be extracted.}
\item{assay}{A string specifying whether counts or normalized expression values are to be extracted.}
\item{value}{
A character vector specifying which control sets are spike-ins.
Alternatively a \code{NULL} value, to remove existing spike-in specifications.
}
}

\details{
The \code{isSpike} and \code{spikes} methods indicate which rows correspond to spike-ins and their expression values, respectively.
If multiple spike-in sets are available, users can extract information for specific sets by supplying the names of the set in \code{type}.
(By default, rows from all spike-in sets are extracted when \code{type=NULL}.)
If \code{assay="exprs"}, users should have run \code{x} through \code{\link{normalize}}.

To specify rows as corresponding to spike-ins, we assume that \code{\link[scater]{calculateQCMetrics}} has already been applied to \code{x}.
Specifically, we assume that spike-ins represent a subset of the control sets supplied as \code{feature_controls} in \code{\link[scater]{calculateQCMetrics}}.
Users can assign a character vector to \code{isSpike(x)<-} containing the names of the control sets that are spike-ins.
This will automatically construct a logical vector containing rows from all specified sets, for later retrieval with \code{isSpike(x)}.

Setting \code{isSpike(x)<-NULL} will clear all existing spike-in information in \code{x}.
Note that direct assignment of a logical vector to \code{isSpike(x)<-} is no longer permitted.
This is because the names of the spike-in sets are necessary for downstream processing, but will not be included if \code{isSpike(x)} is set directly.
}

\value{
For \code{spikes}, a numeric matrix of counts or normalized expression values, with one column per cell and one row per spike-in transcript.

For \code{isSpike}, a logical vector indicating which rows are spike-ins (or \code{NULL}, if this information was not stored in \code{x}).

For \code{isSpike<-}, \code{x} is modified to store a spike-specifying vector in \code{fData(x)$is_feature_spike}.
A logical vector indicating which controls are spike-ins is also stored in the \code{featureControlInfo} slot of \code{x}.
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
y <- calculateQCMetrics(y, list(IAmASpike=rep(c(FALSE, TRUE), c(ngenes, 100))))
isSpike(y) <- "IAmASpike"

y <- computeSpikeFactors(y)
y <- normalize(y)
spikes(y)[1:10,]
spikes(y, assay="exprs")[1:10,]
isSpike(y)
}
