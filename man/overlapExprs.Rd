\name{overlapExprs}
\alias{overlapExprs}
\alias{overlapExprs,matrix-method}
\alias{overlapExprs,SCESet-method}

\title{Overlap expression profiles}
\description{Compute the gene-specific overlap in expression profiles between two groups of cells.}

\usage{
\S4method{overlapExprs}{matrix}(x, groups, design=NULL, residuals=FALSE, tol=1e-8, subset.row=NULL)

\S4method{overlapExprs}{SCESet}(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix of expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SCESet object containing such a matrix.
}
\item{groups}{
A vector of group assignments for all cells.
}
\item{design}{A numeric matrix containing blocking terms, i.e., uninteresting factors driving expression across cells.}
\item{residuals}{A logical scalar indicating whether overlaps should be computed between residuals of a linear model.} 
\item{tol}{A numeric scalar specifying the tolerance with which ties are considered.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to the matrix method.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\details{
For two groups of cells A and B, consider the distribution of expression values for gene X across those cells.
The overlap proportion is defined as the probability that a randomly selected cell in A has a greater expression value of X than a randomly selected cell in B.
Overlap proportions near 0 or 1 indicate that the expression distributions are well-separated.
In particular, large proportions indicate that most cells of the first group (A) express the gene more highly than most cells of the second group (B).

This function computes, for each gene, the overlap proportions between all pairs of groups in \code{groups}.
It is designed to complement \code{\link{findMarkers}}, which reports the log-fold changes between groups.
This is useful for prioritizing candidate markers that are distinctive to one group or another, without needing to plot the expression values.

Expression values that are tied between groups are considered to be 50\% likely to be greater in either group.
Thus, if all values were tied, the overlap proportion would be 0.5.
The tolerance with which ties are considered can be set by changing \code{tol}.

By default, spike-in transcripts are ignored in \code{overlapExprs,SCESet-method} with \code{get.spikes=FALSE}.
This is overridden by any non-\code{NULL} value of \code{subset.row}.
}

\section{Accounting for uninteresting variation}{
If the experiment has known (and uninteresting) factors of variation, these can be included in \code{design}.
The approach used to remove these factors depends on the design matrix.
If there is only one factor in \code{design}, the levels of the factor are defined as separate blocks. 
Overlaps between groups are computed within each block, and a weighted mean (based on the number of cells in each block) of the overlaps is taken across all blocks.

This approach avoids the need for linear modelling and the associated assumptions regarding normality and correct model specification.
In particular, it avoids problems with breaking of ties when counts or expression values are converted to residuals.
However, it also makes less use of information, e.g., we ignore any blocks containing cells from only one group.
\code{NA} proportions may be observed for a pair of groups if there is no block that contains cells from that pair.

For designs containing multiple factors or covariates, a linear model is fitted to the (log-normalized) expression values with \code{design}.
The overlaps are then computed from the residuals of the fitted model.
This approach is not ideal (see \code{?\link{correlatePairs}} for a related discussion) but is unavoidable for covariates.
It can also be used for one-way layouts by setting \code{residuals=TRUE}.
}

\value{
A named list of numeric matrices.
Each matrix corresponds to a group (A) in \code{groups} and contains one row per gene in \code{x} (or the subset specified by \code{subset.row}).
Each column corresponds to another group (B) in \code{groups}.
The matrix entries contain overlap proportions between groups A and B for each gene.
} 

\author{
Aaron Lun
}

\seealso{
\code{\link{findMarkers}}
}

\examples{
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.
groups <- sample(3, ncol(y), replace=TRUE)
out <- overlapExprs(y, groups, subset.row=1:10)
}
