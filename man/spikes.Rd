\name{Get spikes}
\alias{spikes}
\alias{spikes,SCESet-method}

\alias{isSpike}
\alias{isSpike,SCESet-method}
\alias{whichSpike}
\alias{whichSpike,SCESet-method}

\alias{isSpike<-}
\alias{isSpike<-,SCESet-method}
\alias{setSpike<-}
\alias{setSpike<-,SCESet,character-method}
\alias{setSpike<-,SCESet,NULL-method}


\title{Construct the spike-in matrix}
\description{Identify rows in the SCESet corresponding to spike-in transcripts, and retrieve a matrix of counts or normalized expression values for those rows.}

\usage{
\S4method{isSpike}{SCESet}(x, type=NULL)
\S4method{whichSpike}{SCESet}(x)
\S4method{spikes}{SCESet}(x, assay="counts", type=NULL)

\S4method{setSpike}{SCESet,character}(x) <- value
\S4method{setSpike}{SCESet,NULL}(x) <- value
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
If multiple spike-in sets are available, users can extract information for specific sets by supplying the names of the desired set in \code{type}.
(By default, rows from all spike-in sets are extracted when \code{type=NULL}.)
If \code{assay="exprs"}, users should have run \code{x} through \code{\link{normalize}}.
The names of the spike-in sets can also be retrieved with \code{whichSpike}.

To specify rows as corresponding to spike-ins, we assume that \code{\link[scater]{calculateQCMetrics}} has already been applied to \code{x}.
Specifically, we assume that spike-ins represent a subset of the control sets supplied as \code{feature_controls} in \code{\link[scater]{calculateQCMetrics}}.
Users can assign a character vector to \code{setSpike(x)<-} containing the names of the control sets that are spike-ins.
This will automatically construct a logical vector containing rows from all specified sets, for later retrieval with \code{isSpike(x)}.

Setting \code{setSpike(x)<-NULL} will clear all existing spike-in information in \code{x}.
Note that direct assignment of a logical vector to \code{setSpike(x)<-} is not permitted.
This is because the names of the spike-in sets are necessary for downstream processing, but will not be included if \code{isSpike(x)} is set directly.
}

\section{Note on overlapping sets}{
While it is possible to declare overlapping sets as the spike-in sets with \code{isSpike(x)<-}, this is not advisable.
This is because some downstream operations assume that each row belongs to only one set (i.e., one of the spike-in sets, or the set of endogenous genes).
For example, normalization will use size factors from only one of the sets, so correspondence to multiple sets will not be honoured.
A warning will thus be raised if overlapping sets are specified in \code{value}.
}

\value{
For \code{spikes}, a numeric matrix of counts or normalized expression values, with one column per cell and one row per spike-in transcript.

For \code{isSpike}, a logical vector indicating which rows are spike-ins (or \code{NULL}, if this information was not stored in \code{x}).

For \code{whichSpike}, a character vector containing the names of the spike-in sets.

For \code{setSpike<-}, \code{x} is modified to store a spike-specifying vector in \code{fData(x)$is_feature_spike}.
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
example(computeSpikeFactors) # Using the mocked-up data 'y' in this example

# Setting the spike-ins and computing normalized values for them.
setSpike(y) <- "Spike"
y <- computeSpikeFactors(y) # can also do with general.use=FALSE.
y <- normalize(y) 

# Extracting the spike-in data.
head(spikes(y))
head(spikes(y, assay="exprs"))
whichSpike(y)
summary(isSpike(y))
}
