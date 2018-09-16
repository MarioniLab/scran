\name{combinePValues}
\alias{combinePValues}

\title{Combine p-values}
\description{Combine p-values from independent or dependent hypothesis tests using a variety of meta-analysis methods.}

\usage{
combinePValues(..., method=c("fisher", "z", "simes", "berger"),
    weights=NULL, log.p=FALSE)
}

\arguments{
\item{...}{Two or more numeric vectors of p-values of the same length.}
\item{method}{A string specifying the combining strategy to use.}
\item{weights}{A numeric vector of positive weights, with one value per vector in \code{...}.
Alternatively, a list of numeric vectors of weights, with one vector per element in \code{...}.
This is only used when \code{method="z"}.}
\item{log.p}{Logical scalar indicating whether the p-values in \code{...} are log-transformed.}
}

\details{
This function will operate across elements on \code{...} in parallel to combine p-values.
That is, the set of first p-values from all vectors will be combined, followed by the second p-values and so on.
This is useful for combining p-values for each gene across different hypothesis tests.

Fisher's method, Stouffer's Z method and Simes' method test the joint null hypothesis, i.e., that all of the individual null hypotheses in the set are true.
\itemize{
    \item Fisher's method requires independence of the test statistic.
        It is useful in asymmetric scenarios, i.e., when the null is only rejected in one of the tests in the set.
    \item Stouffer's Z method require independence of the test statistic.
        It favours symmetric rejection and can accommodate weights.
    \item Simes' method technically requires independence but tends to be quite robust to dependencies between tests.
}

Berger's intersection-union test examines a different global null hypothesis, i.e., that at least one of the individual null hypotheses are true.
Rejection in the IUT indicates that all of the individual nulls have been rejected.
This is a statistically rigorous equivalent of a naive intersection operation.
However, note that it can be quite conservative.
}

\value{
A numeric vector containing the combined p-values.
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
}

\examples{
p1 <- runif(10000)
p2 <- runif(10000)
p3 <- runif(10000)

fish <- combinePValues(p1, p2, p3)
hist(fish)

z <- combinePValues(p1, p2, p3, method="z", weights=1:3)
hist(z)

simes <- combinePValues(p1, p2, p3, method="simes")
hist(simes)

berger <- combinePValues(p1, p2, p3, method="berger")
hist(berger)
}
