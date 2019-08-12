\name{multiBlockNorm}
\alias{multiBlockNorm}

\title{Per-block scaling normalization}
\description{Perform scaling normalization within each block in a manner that preserves the relative scale of counts between spike-in transcripts and endogenous genes.

This has now been deprecated as rescaling is performed automatically in \code{\link{modelGeneVarWithSpikes}}.
}

\usage{
multiBlockNorm(x, block, ...) 
}

\arguments{
\item{x}{A SingleCellExperiment object containing counts and size factors.}
\item{block}{A factor specifying the blocking level for each cell in \code{x}.}
\item{...}{Further arguments to pass to \code{\link[scater]{normalize}}.}
}

\details{
When comparing spike-in and endogenous variances, the assumption is that a spike-in transcript is comparable to an endogenous gene with similar magnitudes of the counts.
This is motivated by the mean-variance relationship in Poisson and other count-based models.
Thus, we want to ensure that the (relative) average abundances after normalization reflect the similarity in the average count between spike-in transcripts and endogenous genes.
This is usually achieved by centering all sets of size factors so that the normalization does not systematically alter the mean between spike-ins and endogenous genes.
Indeed, this is the default mode of operation in \code{\link[scater]{normalize}}.

However, centering across all cells is not appropriate when \code{block} contains multiple levels and we want to fit trends within each level (see \code{\link{multiBlockVar}}).
In such cases, we want size factors to be centered \emph{within} each level, which is not guaranteed by global centering.
To overcome this, we adjust the spike-in size factors so that the mean within each level of the blocking factor is the same as that of the endogenous size factors for that level.
This avoids cases where spike-in abundances are systematically shifted up or down relative to the abundances of the endogenous genes 
(e.g., due to addition of different spike-in quantities across blocks).

In all cases, the outcome of normalization for endogenous genes is guaranteed to be the same as that from \code{\link[scater]{normalize}}.
Only the size factors and normalized values for spike-in transcripts will be different when using this function.
This ensures that comparisons of gene-level expression profiles between cells (e.g., during clustering or dimensionality reduction) are not altered.
}

\value{
A SingleCellExperiment with normalized log-expression values in the \code{"logcounts"} slot (depending on the arguments to \code{\link[scater]{normalize}}.    
}

\author{
Aaron Lun
}

\seealso{
\code{\link[scater]{normalize}},
\code{\link{multiBlockVar}}
}

\examples{
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.

# Normalizing (gene-based factors for genes, spike-in factors for spike-ins)
y <- computeSumFactors(y) 
y <- computeSpikeFactors(y, general.use=FALSE)

# Setting up the blocking levels.
block <- sample(3, ncol(y), replace=TRUE)
y <- multiBlockNorm(y, block)
assayNames(y)
}

\keyword{normalization}
