\name{testVar}
\alias{testVar}

\title{Test for significantly large variances}
\description{Test for whether the total variance exceeds that expected under some null hypothesis, for sample variances estimated from normally distributed observations.}

\usage{
testVar(total, null, df, design=NULL, test=c("chisq", "f"), 
    fit=NULL, verbose=FALSE)
}

\arguments{
\item{total}{
A numeric vector of total variances for all genes.
}
\item{null}{
A numeric scalar or vector of expected variances under the null hypothesis for all genes.
}
\item{df}{
An integer scalar specifying the degrees of freedom on which the variances were estimated.
}
\item{design}{
A design matrix, used to determine the degrees of freedom if \code{df} is missing.
}
\item{test}{
A string specifying the type of test to perform.
}
\item{fit}{
A list object produced by \code{\link{trendVar}}, only required if \code{test="f"}.
}
\item{verbose}{
A logical scalar specifying whether test parameters should be reported.
}
}

\details{
The null hypothesis is that the true variance for each gene is equal to \code{null}.
(Technically, it is that the variance is equal to or less than this value, but the most conservative test is obtained at equality.)
If \code{test="chisq"}, variance estimates are assumed to follow a chi-squared distribution on \code{df} degrees of freedom and scaled by \code{null/df}.
This is used to compute a p-value for \code{total} being greater than \code{null}.
The underlying assumption is that the observations are normally distributed under the null, which is reasonable for log-counts with low-to-moderate dispersions.

% Also protects against outliers, whereas using the CV2 (e.g., in Brennecke's test method) doesn't.

The aim is to use this function to identify significantly highly variable genes (HVGs).
For example, the \code{null} vector can be set to the values of the trend fitted to the spike-in variances.
This will identify genes with variances significantly greater than technical noise.
Alternatively, it can be set to the trend fitted to the cellular variances, which will identify those that are significantly more variable than the bulk of genes.
Selecting HVGs on p-values is better than using \code{total - null}, as the latter is less precise when \code{null} is large.

If \code{test="f"}, the true variance of each spike-in transcript is assumed to be sampled from a scaled inverse chi-squared distribution.
This accounts for any inflated scatter around the trend due to differences in amplification efficiency between transcripts.
As a result, the gene-wise variance estimates are F-distributed under the null, which is exploited to calculate the p-value for each gene.
Parameters of the F distribution are estimated from values in \code{fit} using \code{\link[limma]{fitFDistRobustly}}.

% Each transcript is treated as a random variable; thus, while error may not be controlled for a specific transcript, it will be controlled over all transcripts.
}

\value{
A numeric vector of p-values for all genes.
}

\seealso{
\code{\link{trendVar}},
\code{\link{decomposeVar}},
\code{\link[limma]{fitFDistRobustly}}
}

\author{
Aaron Lun
}

\examples{
set.seed(100)
null <- 100/runif(1000, 50, 2000)
df <- 30
total <- null * rchisq(length(null), df=df)/df

# Direct test:
out <- testVar(total, null, df=df)
hist(out)

# Rejecting the null:
alt <- null * 5 * rchisq(length(null), df=df)/df
out <- testVar(alt, null, df=df)
plot(alt[order(out)]-null)

# Focusing on genes that have high absolute increases in variability:
out <- testVar(alt, null+0.5, df=df)
plot(alt[order(out)]-null)
}

\keyword{variance}

\references{
Law CW, Chen Y, Shi W and Smyth GK (2014).
voom: precision weights unlock linear model analysis tools for RNA-seq read counts
\emph{Genome Biol.} 15(2), R29. 
}

