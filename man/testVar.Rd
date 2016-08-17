\name{testVar}
\alias{testVar}

\title{Test for significantly large variances}
\description{Test for whether the total variance exceeds that expected under some null hypothesis, for sample variances estimated from normally distributed observations.}

\usage{
testVar(total, null, df, design=NULL) 
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
}

\details{
The null hypothesis states that the true variance for each gene is equal to \code{null}.
(Technically, it states that the variance is equal to or less than this value, but the most conservative test is obtained at equality.)
Variance estimates should follow a scaled chi-squared distribution on \code{df} degrees of freedom, where the scaling factor is equal to \code{null/df}.
This can be used to compute a p-value for \code{total} being greater than \code{null}.
The assumption is that the original observations were normally distributed -- using log-CPMs tends to work reasonably well for count data.

% Also protects against outliers, whereas using the CV2 (e.g., in Brennecke's test method) doesn't.

The idea is to use this function to identify significantly highly variable genes (HVGs).
For example, the \code{null} vector can be set to the values of the trend fitted to the spike-in variances.
This will identify genes with variances significantly greater than technical noise.
Alternatively, it can be set to the trend fitted to the cellular variances, which will identify those that are significantly more variable than the bulk of genes.
Ranking HVGs on p-values is better than ranking on \code{total - null}, as the latter is less precise when \code{null} is large.
}

\value{
A numeric vector of p-values for all genes.
}

\seealso{
\code{\link{trendVar}},
\code{\link{decomposeVar}}
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

