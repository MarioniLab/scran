\name{combineVar}
\alias{combineVar}

\title{Combine variance decompositions}
\description{Combine the results of multiple variance decompositions, usually generated for the same genes across separate batches of cells.}

\usage{
combineVar(..., method="fisher", weighted=TRUE)
}

\arguments{
\item{...}{Two or more DataFrames produced by \code{\link{decomposeVar}}.}
\item{method}{String specifying how p-values are to be combined, see \code{\link{combinePValues}} for options.}
\item{weighted}{Logical scalar indicating whether weights should be used for combining statistics.}
}

\details{
This function is designed to merge results from multiple calls to \code{\link{decomposeVar}}, usually computed for different batches of cells.
Separate variance decompositions are necessary in cases where different concentrations of spike-in have been added to the cells in each batch.
This affects the technical mean-variance relationship and precludes the use of a common trend fit.

The output mean is computed as a weighted average of the means in each input DataFrame, where the weight is defined as the number of cells in that batch.
This yields an equivalent value to the sample mean across all cells in all batches.
Similarly, weighted averages are computed for all variance components, where the weight is defined as the residual d.f. used for variance estimation in each batch.
This yields a variance equivalent to the residual variance obtained while blocking on the batch of origin.

Weighting can be turned off with \code{weighted=FALSE}.
This may be useful to ensure that all batches contribute equally to the calculation of the combined statistics, avoiding cases where batches with many cells dominate the output.
Of course, this comes at the cost of precision - large batches contain more information and \emph{should} contribute more to the weighted average.

For the p-value calculations, options are taken from \code{\link{combinePValues}} and are paraphrased here:
\itemize{
    \item The default approach is to use \code{method="fisher"}, which will identify genes that are highly variable in \emph{any} batch.
    \item Another option is to use \code{method="z"}, which favours (but does not require) genes that are significant in multiple batches.
    \item Setting \code{method="simes"} is a more conservative counterpart to Fisher's method that is robust to correlations.
    \item Setting \code{method="berger"} will identify genes that are significant in \emph{all} batches.
}

Only \code{method="z"} will perform any weighting of batches, and only if \code{weighted=TRUE}.
Here, each batch is weighted according to the residual d.f. used for testing in that batch.
In all other cases, all batches are assigned equal weight for p-value calculations.
}

\value{
A DataFrame with the same numeric fields as that produced by \code{\link{decomposeVar}}.
Each field contains the average across all batches except for \code{p.value}, which contains the combined p-value based on \code{method};
and \code{FDR}, which contains the adjusted p-value using the BH method.
}

\seealso{
\code{\link{decomposeVar}},
\code{\link{combinePValues}}
}

\author{
Aaron Lun
}

\examples{
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.
y <- computeSumFactors(y) # Size factors for the the endogenous genes.
y <- computeSpikeFactors(y, general.use=FALSE) # Size factors for spike-ins. 

y1 <- y[,1:100] 
y1 <- normalize(y1) # normalize separately after subsetting.
fit1 <- trendVar(y1)
results1 <- decomposeVar(y1, fit1)

y2 <- y[,1:100 + 100] 
y2 <- normalize(y2) # normalize separately after subsetting.
fit2 <- trendVar(y2)
results2 <- decomposeVar(y2, fit2)

head(combineVar(results1, results2))
head(combineVar(results1, results2, method="simes"))
head(combineVar(results1, results2, method="berger"))
}
