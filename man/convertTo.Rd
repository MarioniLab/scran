\name{convertTo}
\alias{convertTo}
\alias{convertTo,SCESet-method}

\title{Convert to other classes}
\description{
Convert a SCESet object into other classes for entry into other analysis pipelines.
}

\usage{
\S4method{convertTo}{SCESet}(x, type=c("edgeR", "DESeq2", "monocle"),
    fData.col=NULL, pData.col=NULL, ..., assay, 
    use.all.sf=TRUE, normalize=TRUE, subset.row=NULL, get.spikes=FALSE)
}

\arguments{
\item{x}{A SCESet object.}
\item{type}{A string specifying the analysis for which the object should be prepared.}
\item{fData.col}{Any set of indices specifying which columns of \code{fData(x)} should be retained in the returned object.}
\item{pData.col}{Any set of indices specifying which columns of \code{pData(x)} should be retained.}
\item{...}{Other arguments to be passed to pipeline-specific constructors.}
\item{assay}{A string specifying which assay of \code{x} should be put in the returned object.}
\item{use.all.sf}{A logical scalar indicating whether multiple size factors should be used to generate the returned object.}
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
This is consistent with the required inputs to these analyses.
Information about normalization is instead transmitted via size or normalization factors in the output object.
For monocle, \code{assay} is ignored and counts are divided by the size factors to yield (roughly) log-normally distributed expression values.
Values in \code{assay} can be used directly by setting \code{normalize=FALSE}.

In all cases, rows corresponding to spike-in transcripts are removed from the output object by default.
As such, rows in the returned object may not correspond directly to rows in \code{x}.
Users should consider this when retrieving analysis results from these pipelines, e.g., match on row names in \code{x} before comparing to other results.
This behaviour can be turned off by setting \code{get.spikes=TRUE}, such that all rows are retrieved in the output object.
Users can also set \code{subset.row} to extract specific rows, in which case \code{get.spikes} is ignored.

By default, different size factors for different rows (e.g., for spike-in sets) will be respected.
For edgeR, an offset matrix will be constructed containing mean-centred log-size factors for each row.
For DESeq2, a similar matrix will be constructed containing size factors scaled to have a geometric mean of unity.
For monocle, counts for each row will be divided by the size factors for that row.
This behaviour can be turned off with \code{use.all.sf=FALSE}, such that only \code{sizeFactors(x)} is used for normalization for all \code{type}.
(For edgeR and DESeq2, the offset matrix is not generated if all rows correspond to \code{sizeFactors(x)}, as this information is already stored in the object.)
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
example(computeSpikeFactors) # Using the mocked up data 'y' from this example.
sizeFactors(y) <- 2^rnorm(ncells) # Adding some additional embellishments.
fData(y)$SYMBOL <- paste0("X", seq_len(nrow(y)))
y$other <- sample(LETTERS, ncells, replace=TRUE)

# Converting to various objects.
convertTo(y, type="edgeR")
convertTo(y, type="DESeq2")
convertTo(y, type="monocle")
}

