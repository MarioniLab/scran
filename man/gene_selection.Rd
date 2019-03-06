\name{Gene selection}
\alias{scran-gene-selection}

\title{Gene selection}
\description{Details on how gene selection is performed in almost all \pkg{scran} functions.}

\section{Subsetting by row}{ 
For functions accepting some matrix \code{x}, we can choose to perform calculations only on a subset of rows with \code{subset.row}.
This can be a logical, integer or character vector indicating the rows of \code{x} to use.
If a character vector, it must contain the names of the rows in \code{x}.
Future support will be added for more esoteric subsetting vectors like the Bioconductor \linkS4class{Rle} classes.

The output of running a function with \code{subset.row} will \emph{always} be the same as the output of subsetting \code{x} beforehand and passing it into the function.
However, it is often more efficient to use \code{subset.row} as we can avoid constructing an intermediate subsetted matrix.
The same reasoning applies for any \code{x} that is a \linkS4class{SingleCellExperiment} object.
}

\section{Handling spike-in transcripts}{
Many functions only make sense when performed on the endogenous genes.
For such functions, spike-in transcripts are automatically removed when the input is a SingleCellExperiment and \code{get.spikes=FALSE}.
This is achieved based on the spike-in information returned by \code{\link{isSpike}}.

If \code{get.spikes=TRUE}, no filtering on the spike-in transcripts will be performed.
Filtering will not be performed (and in fact, the argument should be unavailable) if the input is not a SingleCellExperiment object.

If \code{get.spikes=FALSE} and \code{subset.row} is not \code{NULL}, the two selections are intersected.
That is, only the non-spike-in entries of \code{subset.row} will be used in the function.

A number of functions may require special treatment of spike-ins, e.g., \code{\link{trendVar}}.
Refer to the corresponding documentation for more details.
}

\section{Filtering by mean}{
Some functions will have a \code{min.mean} argument to filter out low-abundance genes prior to processing.
Depending on the function, the filter may be applied to the average library size-adjusted count computed by \code{\link{calcAverage}}, the average log-count, 
or some other measure of abundance - see the documentation for each function for details.

Any filtering on \code{min.mean} is automatically intersected with \code{get.spikes=FALSE} and/or a specified \code{subset.row}.
For example, only non-spike-in genes that pass the \code{min.mean} filter are retained if \code{get.spikes=TRUE}.
Similarly, only selected genes that pass the filter are retained if \code{subset.row} is specified.
}

\author{
Aaron Lun
}
