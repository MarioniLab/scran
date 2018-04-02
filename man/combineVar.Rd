\name{combineVar}
\alias{combineVar}

\title{Combine variance decompositions}
\description{Combine the results of multiple variance decompositions, usually generated for the same genes across separate batches of cells.}

\usage{
combineVar(..., method=c("fisher", "z", "simes", "berger"), weighted=TRUE)
}

\arguments{
\item{...}{Two or more DataFrames produced by \code{\link{decomposeVar}}.}
\item{method}{String specifying how p-values are to be combined.}
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
}

\section{Combining p-values across batches}{
The default approach is to use \code{method="fisher"}, which uses Fisher's method for combining p-values.
This will test the global null hypothesis, i.e., that genes are not variable in any batch.
As a result, it will identify genes that are highly variable in \emph{any} batch.

Another option is to use \code{method="z"}, where Stouffer's Z-score method is used to combine p-values across batches.
Each batch is weighted according to the residual d.f. used to perform the test.
This is less sensitive that Fisher's method to low p-values in only a single batch, instead favouring genes that are significant in multiple batches.

Both Fisher's and Stouffer's methods assume independence between batches. 
If this is not the case, Simes' method should be used by setting \code{method="simes"}.
This is more robust to correlations between tests (Sarkar and Chung, 1997; see also similar work on the related Benjamini-Hochberg method).

To identify genes that are detected as highly variable in \emph{all} batches, Berger's IUT can be used by setting \code{method="berger"}.
This defines the combined p-value as the maximum across batches for each gene.
Needless to say, this is quite a conservative approach.

Only \code{method="z"} will perform any weighting of batches, and only if \code{weighted=TRUE}.
In all other cases, all batches are assigned equal weight.
}

\value{
A DataFrame with the same numeric fields as that produced by \code{\link{decomposeVar}}.
Each field contains the average across all batches except for \code{p.value}, which contains the combined p-value based on \code{method};
and \code{FDR}, which contains the adjusted p-value using the BH method.
}

\seealso{
\code{\link{decomposeVar}}
}

\author{
Aaron Lun
}

\references{
Fisher, R.A. (1925). 
\emph{Statistical Methods for Research Workers.}
Oliver and Boyd (Edinburgh).

Whitlock MC (2005). 
Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
\emph{J. Evol. Biol.} 18, 5:1368-73.

Simes RJ (1986). 
An improved Bonferroni procedure for multiple tests of significance. 
\emph{Biometrika} 73:751-754.

Berger RL and Hsu JC (1996). 
Bioequivalence trials, intersection-union tests and equivalence confidence sets.
\emph{Statist. Sci.} 11, 283-319.

Sarkar SK and Chung CK (1997).
The Simes method for multiple hypothesis testing with positively dependent test statistics.
\emph{J. Am. Stat. Assoc.} 92, 1601-1608.
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
