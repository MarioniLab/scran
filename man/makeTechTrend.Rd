\name{makeTechTrend}
\alias{makeTechTrend}

\title{Make a technical trend}
\description{Manufacture a mean-variance trend for log-transformed expression values, assuming Poisson or NB-distributed technical noise for count data.}

\usage{
makeTechTrend(means, size.factors=1, tol=1e-6, dispersion=0, 
    pseudo.count=1, sce=NULL) 
}

\arguments{
\item{means}{A numeric vector of average counts.
Note that there are means of the counts, \emph{not} means of the log-expression values.}
\item{size.factors}{A numeric vector of size factors.
If supplied, these should be centred at unity.}
\item{tol}{A numeric scalar specifying the tolerance for approximating the mean/variance.
Lower values result in greater accuracy.}
\item{dispersion}{A numeric scalar specifying the dispersion for the NB distribution.
If zero, this is equivalent to a Poisson distribution.}
\item{pseudo.count}{A numeric scalar specifying the pseudo-count to be added to the scaled counts before log-transformation.}
\item{sce}{A SingleCellExperiment object from which \code{size.factors} and \code{pseudo.count} are extracted,
and \code{means} can be automatically inferred.}
}

\details{
At each value of \code{means}, this function will examine the distribution of Poisson/NB-distributed counts with the corresponding mean.
All counts are log2-transformed after addition of \code{pseudo.count}, and the mean and variance is computed for the log-transformed values.
Setting \code{dispersion} to a non-zero value will use a NB distribution instead of the default Poisson.

If \code{size.factors} is a vector, one count distribution is generated for each of its elements, where the mean is scaled by the corresponding size factor.
Counts are then divided by the size factor prior to log-transformation, mimicking the effect of normalization in \code{\link[scater]{normalize}}.
A composite distribution of log-values is constructed by pooling all of these individual distributions.
The mean and variance is then computed for a composite distribution.

Finally, a function is fitted to all of the computed variances, using the means of the log-values as the covariate.
Note that the returned function accepts mean log-values as input, \emph{not} the mean counts that were supplied in \code{means}.
This means that the function is directly usable as a replacement for the trend fitted by \code{\link{trendVar}}.

If \code{sce} is set, \code{pseudo.count} is overridden by \code{metadata(sce)$log.exprs.offset};
\code{size.factors} is overridden by \code{sizeFactors(sce)} (or the column sums of \code{counts(sce)}, if no size factors are specified);
and \code{means} is automatically determined from the range of row averages of \code{logcounts(sce)} (after undoing the log-transformation).
}

\value{
A function accepting a mean log-expression as input and returning the variance of the log-expression as the output.
}

\seealso{
\code{\link{trendVar}},
\code{\link[scater]{normalize}}
}

\author{
Aaron Lun
}

\examples{
means <- 1:100/10
out <- makeTechTrend(means)
curve(out(x), xlim=c(0, 5))
}
