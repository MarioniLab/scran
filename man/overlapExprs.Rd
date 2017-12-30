\name{overlapExprs}
\alias{overlapExprs}
\alias{overlapExprs,ANY-method}
\alias{overlapExprs,SingleCellExperiment-method}

\title{Overlap expression profiles}
\description{Compute the gene-specific overlap in expression profiles between two groups of cells.}

\usage{
\S4method{overlapExprs}{ANY}(x, groups, block=NULL, design=NULL, 
    rank.type=c("any", "all"), direction=c("any", "up", "down"),
    tol=1e-8, BPPARAM=SerialParam(), subset.row=NULL, 
    lower.bound=NULL, residuals=FALSE)

\S4method{overlapExprs}{SingleCellExperiment}(x, ..., subset.row=NULL, lower.bound=NULL, 
    assay.type="logcounts", get.spikes=FALSE) 
}

\arguments{
\item{x}{
    A numeric matrix of expression values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{groups}{
A vector of group assignments for all cells.
}
\item{block}{A factor specifying the blocking level for each cell.}
\item{design}{A numeric matrix containing blocking terms, i.e., uninteresting factors driving expression across cells.}
\item{rank.type}{A string specifying which comparisons should be used to rank genes in the output.}
\item{direction}{A string specifying which direction of change in expression should be used to rank genes in the output.}
\item{tol}{A numeric scalar specifying the tolerance with which ties are considered.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{lower.bound}{A numeric scalar specifying the theoretical lower bound of values in \code{x}, only used when \code{residuals=TRUE}.}
\item{residuals}{A logical scalar indicating whether overlaps should be computed between residuals of a linear model.} 
\item{...}{Additional arguments to pass to the matrix method.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\details{
For two groups of cells A and B, consider the distribution of expression values for gene X across those cells.
The overlap proportion is defined as the probability that a randomly selected cell in A has a greater expression value of X than a randomly selected cell in B.
Overlap proportions near 0 or 1 indicate that the expression distributions are well-separated.
In particular, large proportions indicate that most cells of the first group (A) express the gene more highly than most cells of the second group (B).

This function computes, for each gene, the overlap proportions between all pairs of groups in \code{groups}.
It will then rank the genes based on how well they differentiate between groups.
\code{overlapExprs} is designed to complement \code{\link{findMarkers}}, which reports the log-fold changes between groups.
This is useful for prioritizing candidate markers without needing to plot their expression values.

Expression values that are tied between groups are considered to be 50\% likely to be greater in either group.
Thus, if all values were tied, the overlap proportion would be equal to 0.5.
The tolerance with which ties are considered can be set by changing \code{tol}.

Users can specify which subset of genes to perform these calculations on, by supplying a non-\code{NULL} value of \code{subset.row}.
By default, spike-in transcripts are ignored in \code{overlapExprs,SingleCellExperiment-method} with \code{get.spikes=FALSE}.
If \code{get.spikes=FALSE} and \code{subset.row!=NULL}, the function will only use the non-spike-in transcripts in \code{subset.row}.
}

\section{Ranking genes in the output}{
If \code{rank.type="any"}, the genes in each group-specific DataFrame are ranked using a similar logic to that in \code{\link{findMarkers}}.
In the DataFrame for group A, a \code{Top} value is calculated for each gene, representing an aggregated ranking across all pairwise comparisons to the other groups.
To illustrate, consider the set of all genes with \code{Top} values less than or equal to some integer X.
This set is the union of the top X genes with the worst overlaps from each pairwise comparison between group A and every other group.
Ranking genes based on the \code{Top} value prioritizes genes that exhibit low overlaps between group A and \emph{any} other group.

If \code{rank.type="all"}, the genes in each group-specific DataFrame are ranked by the \code{Best} value instead.
This is the proportion corresponding to the best overlap across all pairwise comparisons between group A and the other groups.
By using this metric, genes can only achieve a high ranking if there exist poor overlaps between group A and \emph{all} other groups.
This tends to be quite conservative but can be helpful for quickly identifying uniquely differentially expressed markers.

By default, \code{direction="any"} means that any shift in the distribution between groups that reduces overlap will be considered during ranking.
Thus, poor overlaps are defined by overlap proportions with the largest absolute difference from 0.5.
Conversely, good overlaps are defined by overlap proportions closest to 0.5.
If \code{direction="up"}, the function will instead focus on genes that shift upwards in A versus all other groups, 
i.e., poor overlaps have large overlap proportions and good overlaps have small proportions.
The opposite is true if \code{direction="down"}, where the function will focus on genes that shift downwards in A versus all other groups.
}

\section{Accounting for uninteresting variation}{
If the experiment has known (and uninteresting) factors of variation, these can be included in \code{design} or \code{block}.
The approach used to remove these factors depends on which argument is used.
If there is only one factor, using \code{block} is recommended whereby the levels of the factor are defined as separate groups. 
Overlaps between groups are computed within each block, and a weighted mean (based on the number of cells in each block) of the overlaps is taken across all blocks.

This approach avoids the need for linear modelling and the associated assumptions regarding normality and correct model specification.
In particular, it avoids problems with breaking of ties when counts or expression values are converted to residuals.
However, it also makes less use of information, e.g., we ignore any blocks containing cells from only one group.
\code{NA} proportions may also be observed for a pair of groups if there is no block that contains cells from that pair.

For experiments containing multiple factors or covariates, a linear model is fitted to the expression values with an appropriate matrix in \code{design}.
Overlap proportions are then computed using the residuals of the fitted model.
This approach is not ideal, requiring log-transformed \code{x} and setting of \code{lower.bound} - see \code{?\link{correlatePairs}} for a related discussion.
Where possible for one-way layouts, we suggest using \code{block} instead.
}

\value{
A named list of DataFrames.
Each DataFrame corresponds to a group in \code{groups} and contains one row per gene in \code{x} (or the subset specified by \code{subset.row}).
Within the DataFrame for each group (e.g., group A), there are the following fields:
\describe{
\item{\code{Top}:}{Integer, the aggregated rank across all pairwise comparisons if \code{rank.type="any"}.}
\item{\code{Best}:}{Numeric, the value of the best overlap proportion across all comparisons if \code{rank.type="all"}.}
\item{\code{overlap.B}:}{Numeric for every other group B in \code{groups}, containing overlap proportions between groups A and B for that gene.}
}
Genes are ranked by the \code{Top} or \code{Best} column, depending on \code{rank.type}.
} 

\author{
Aaron Lun
}

\seealso{
\code{\link{findMarkers}}
}

\examples{
# Using the mocked-up data 'y2' from this example.
example(computeSpikeFactors) 
y2 <- normalize(y2)
groups <- sample(3, ncol(y2), replace=TRUE)
out <- overlapExprs(y2, groups, subset.row=1:10)
}
