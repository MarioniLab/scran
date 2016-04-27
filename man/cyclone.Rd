\name{cyclone}
\alias{cyclone}
\alias{cyclone,matrix-method}
\alias{cyclone,SCESet-method}

\title{Cell cycle phase classification}
\description{Classify single cells into their cell cycle phases based on gene expression data.}

\usage{
\S4method{cyclone}{matrix}(x, pairs, gene.names=rownames(x), iter=1000, min.iter=100, min.pairs=50, 
    BPPARAM=bpparam(), verbose=FALSE)
\S4method{cyclone}{SCESet}(x, ..., assay="counts", get.spikes=FALSE)
}

\arguments{
\item{x}{
    A numeric matrix of gene expression values where rows are genes and columns are cells.
    Alternatively, a SCESet object containing such a matrix.
}
\item{pairs}{A list of data.frames produced by \code{\link{sandbag}}, containing pairs of marker genes.}
\item{gene.names}{A character vector of gene names.}
\item{iter}{An integer scalar specifying the number of iterations for random sampling to obtain a cycle score.}
\item{min.iter}{An integer scalar specifying the minimum number of iterations for score estimation.}
\item{min.pairs}{An integer scalar specifying the minimum number of pairs for cycle estimation.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{verbose}{A logical scalar specifying whether diagnostics should be printed to screen.}
\item{...}{Additional arguments to pass to \code{cyclone,matrix-method}.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical specifying whether spike-in transcripts should be used.}
}

\details{
This function implements the classification step of the pair-based prediction method described by Scialdone et al. (2015).
Pairs of marker genes are trained with \code{\link{sandbag}}, where the sign of the relative expression between gene in each pair changes across phases.
For each phase and each cell, the function calculates the proportion of all marker pairs where the expression of the first gene is greater than the second
    (pairs with the same expression are ignored).
A distribution of proportions is constructed by shuffling the expression values within the cell and recalculating the proportion at each iteration.
The phase score for that cell is then defined as the lower tail probability of this distribution.

By default, shuffling is performed \code{iter} times to obtain the distribution from which the score is estimated.
However, some iterations may not be used if there are fewer than \code{min.pairs} pairs with different expression, such that the proportion cannot be calculated precisely.
Also, a score is only returned if the distribution is large enough for stable calculation of the tail probability, i.e., consists of results from at least \code{min.iter} iterations.

% The original implementation shuffled expression values from all genes in the training set, rather than only genes in the marker set for each phase.
% In theory, this shouldn't matter, as the score calculation depends on the ranking of each gene.
% That should be the same regardless of the distribution of expression values -- each set of rankings is equally likely, no matter what.
% In practice, the number of ties will differ between different set of genes.
% The most appropriate comparison would involve the same number of ties as that used to calculate the observed score.
% Thus, I would set the marker set genes as the empirical distribution, rather than the entire training set of genes.

Cells with G1 and G2M scores above 0.5 should be assigned to the G1 and G2M phases, respectively.
This is based on the interpretation of the score as 1 minus the p-value for the null distribution of proportions.
The null hypothesis here is that expression of the marker genes is independent within each cell, i.e., with no cycle-induced correlations between marker pairs.
Cells can be assigned to S phase based on the S phase score, but a more reliable approach is to define S phase cells based on those cells with G1 and G2M scores below 0.5.

% That said, you can have non-cycle-based correlations, so full independence is probably liberal.

For \code{cyclone,SCESet-method}, the matrix of counts is used but can be replaced with expression values by setting \code{assays}.
By default, \code{get.spikes=FALSE} which means that any rows corresponding to spike-in transcripts will not be considered for score calculation.
This is for the same reasons as described in \code{?\link{sandbag}}.
}

\value{
A list of two data frames is returned -- 
    \code{scores}, containing the phase scores for each phase and cell (i.e., each row is a cell); 
    and \code{normalized.scores}, containing the row-normalized scores (i.e., where the row sum for each cell is equal to 1).
}

\author{
Antonio Scialdone,
with modifications by Aaron Lun
}

\seealso{
\code{\link{sandbag}}
}

\examples{
example(sandbag)

# Classifying (note: test.data!=training.data in real cases)
test <- training 
assignments <- cyclone(test, out)

# Visualizing
col <- character(ncells)
col[is.G1] <- "red"
col[is.G2M] <- "blue"
col[is.S] <- "darkgreen"
plot(assignments$score$G1, assignments$score$G2M, col=col, pch=16)
}

\references{
Scialdone A, Natarajana KN, Saraiva LR et al. (2015). 
Computational assignment of cell-cycle stage from single-cell transcriptome data.
\emph{Methods} 85:54--61
}

\keyword{clustering}
