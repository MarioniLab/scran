\name{convertTo}
\alias{convertTo}
\alias{convertTo,SCESet-method}

\title{Convert to other classes}
\description{
Convert a SCESet object into other classes for entry into other analysis pipelines.
}

\usage{
\S4method{convertTo}{SCESet}(x, type=c("edgeR", "DESeq2", "monocle"),
    fData.col=NULL, pData.col=NULL, ..., 
    assay, normalize=TRUE, subset.row=NULL, get.spikes=FALSE)
}

\arguments{
\item{x}{A SCESet object.}
\item{type}{A string specifying the analysis for which the object should be prepared.}
\item{fData.col}{Any set of indices specifying which columns of \code{fData(x)} should be retained in the returned object.}
\item{pData.col}{Any set of indices specifying which columns of \code{pData(x)} should be retained.}
\item{...}{Other arguments to be passed to pipeline-specific constructors.}
\item{assay}{A string specifying which assay of \code{x} should be put in the returned object.}
\item{normalize}{A logical scalar specifying whether the assay values should be normalized for \code{type="monocle"}.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to return.}
\item{get.spikes}{A logical scalar specifying whether rows corresponding to spike-in transcripts should be returned.}
}

\value{
For \code{type="edgeR"}, a DGEList object is returned containing the count matrix.
Size factors are converted to normalization factors.
Gene-specific \code{fData} is stored in the \code{genes} element, and cell-specific \code{pData} is stored in the \code{samples} element.

For \code{type="DESeq2"}, a DESeqDataSet object is returned containing the count matrix and size factors.
Additional gene- and cell-specific data is stored in the \code{mcols} and \code{colData} respectively.

For \code{type="monocle"}, a CellDataSet object is returned containing the unlogged expression values.
Additional gene- and cell-specific data is stored in the \code{fData} and \code{pData} respectively.
}

\details{
This function converts an SCESet object into various other classes in preparation for entry into other analysis pipelines, as specified by \code{type}.
Gene- and cell-specific data fields can be retained in the output object by setting \code{fData.col} and \code{pData.col}, respectively.
Other arguments can be passed to the relevant constructors through the ellipsis.

By default, for edgeR and DESeq2, \code{assay} is set to \code{"counts"} such that count data is stored in the output object.
This is consistent with the required inputs to these analyses, where normalization information is instead transmitted via size or normalization factors.
For monocle, \code{assay} is ignored and counts are divided by the size factors to yield (roughly) log-normally distributed expression values.
Values in \code{assay} can be used directly by setting \code{normalize=FALSE}.

In all cases, rows corresponding to spike-in transcripts are removed from the output object by default.
As such, rows in the returned object may not correspond directly to rows in \code{x}.
Users should consider this when retrieving analysis results from these pipelines, e.g., match on row names in \code{x} before comparing to other results.
This behaviour can be turned off by setting \code{get.spikes=TRUE}, such that all rows are retrieved in the output object.
Users can also set \code{subset.row} to extract specific rows, in which case \code{get.spikes} is ignored.

Note that, for all \code{type}, only the main size factors will be stored or used inside \code{convertTo}.
Any sets of size factors specific to, say, spike-in transcripts (see \code{?\link[scater]{normalize,SCESet-method}}) will be ignored.
}

\author{
Aaron Lun
}

\seealso{
\code{\link[edgeR]{DGEList}},
\code{\link[DESeq2]{DESeqDataSetFromMatrix}},
\code{\link[monocle]{newCellDataSet}}
}

\examples{
ncells <- 200
ngenes <- 100
count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- newSCESet(countData=data.frame(dummy))
is.spike <- rbinom(ngenes, 1, 0.5)==0L
isSpike(X) <- is.spike
sizeFactors(X) <- 2^rnorm(ncells)
X <- normalize(X)

fData(X)$SYMBOL <- paste0("X", seq_len(ngenes))
X$other <- sample(LETTERS, ncells, replace=TRUE)

convertTo(X, type="edgeR")
convertTo(X, type="DESeq2")
convertTo(X, type="monocle")
}

