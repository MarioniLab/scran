\name{Gene selection}
\alias{scran-gene-selection}

\title{Gene selection}
\description{Details on how gene selection is performed in almost all \pkg{scran} functions.}

\section{Subsetting by row}{ 
For functions accepting some gene-by-cell matrix \code{x}, we can choose to perform calculations only on a subset of rows (i.e., genes) with the \code{subset.row} argument.
This can be a logical, integer or character vector indicating the rows of \code{x} to use.
If a character vector, it must contain the names of the rows in \code{x}.
Future support will be added for more esoteric subsetting vectors like the Bioconductor \linkS4class{Rle} classes.

The output of running a function with \code{subset.row} will \emph{always} be the same as the output of subsetting \code{x} beforehand and passing it into the function.
However, it is often more efficient to use \code{subset.row} as we can avoid constructing an intermediate subsetted matrix.
The same reasoning applies for any \code{x} that is a \linkS4class{SingleCellExperiment} object.
}

\section{Filtering by mean}{
Some functions will have a \code{min.mean} argument to filter out low-abundance genes prior to processing.
Depending on the function, the filter may be applied to the average library size-adjusted count computed by \code{\link{calculateAverage}}, the average log-count, 
or some other measure of abundance - see the documentation for each function for details.

Any filtering on \code{min.mean} is automatically intersected with a specified \code{subset.row}.
For example, only subsetted genes that pass the filter are retained if \code{subset.row} is specified alongside \code{min.mean}.
}

\author{
Aaron Lun
}
