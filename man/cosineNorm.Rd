\name{cosineNorm}
\alias{cosineNorm}

\title{Cosine normalize}
\description{Perform cosine-normalization on the column vectors of an expression matrix.}

\usage{
cosineNorm(X, mode=c("matrix", "all", "l2norm"))
}

\arguments{
    \item{X}{A gene expression matrix with cells as columns and genes as rows.}
    \item{mode}{A string specifying the output to be returned.}
}

\details{
While the default is to directly return the cosine-normalized matrix, it may occasionally be desirable to obtain the L2 norm, 
e.g., to apply an equivalent normalization to other matrices.
This can be achieved by setting \code{mode} accordingly.
}

\value{
If \code{mode="matrix"}, a double-precision matrix of the same dimensions as \code{X} is returned, containing cosine-normalized values.

If \code{mode="l2norm"}, a double-precision vector is returned containing the L2 norm for each cell.

If \code{mode="all"}, a named list is returned containing the fields \code{"matrix"} and \code{"l2norm"}, which are as described above.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{mnnCorrect}},
\code{\link{fastMNN}}
}

\examples{
A <- matrix(rnorm(1000), nrow=10)
str(cosineNorm(A))
str(cosineNorm(A, mode="l2norm"))
}

\keyword{normalization}
