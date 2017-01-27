\name{combineVar}
\alias{combineVar}

\title{Combine variance decompositions}
\description{Combine the results of multiple variance decompositions, usually generated for the same genes across separate batches of cells.}

\usage{
combineVar(..., method=c("fisher", "simes", "berger"))
}

\arguments{
\item{...}{Two or more data frames, each produced by \code{\link{decomposeVar}}.}
\item{method}{A string specifying how p-values are to be combined.}
}

\details{
This function is designed to merge results from multiple calls to \code{\link{decomposeVar}}, usually computed for different batches of cells.
Separate variance decompositions are necessary in cases where different concentrations of spike-in have been added to the cells in each batch.
This affects the technical mean-variance relationship and precludes the use of a common trend fit.

The default setting is to use \code{method="fisher"}, where Fisher's method is used to combine p-values across batches.
This aims to detect genes that are highly variable in \emph{any} batch and assumes that the test outcome is independent between batches.
If independence does not hold, Simes' method should be used by setting \code{method="simes"}, as it is more robust to correlations between tests.
To identify genes that are detected as highly variable in \emph{all} batches, Berger's IUT can be used by setting \code{method="simes"}.
}

\value{
A data frame with the same numeric fields as that produced by \code{\link{decomposeVar}}.
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
Simes RJ (1986). 
An improved Bonferroni procedure for multiple tests of significance. 
\emph{Biometrika} 73:751-754.

Berger RL and Hsu JC (1996). 
Bioequivalence trials, intersection-union tests and equivalence confidence sets.
\emph{Statist. Sci.} 11, 283-319.

Fisher, R.A. (1925). 
\emph{Statistical Methods for Research Workers.}
Oliver and Boyd (Edinburgh).
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
