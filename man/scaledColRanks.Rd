\name{Scaled column ranks}
\alias{scaledColRanks}

\title{Compute scaled column ranks}
\description{Compute scaled column ranks from each cell's expression profile for distance calculations based on rank correlations.}

\usage{
scaledColRanks(x, subset.row=NULL, min.mean=NULL, transposed=FALSE, 
    as.sparse=FALSE, withDimnames=TRUE)
}

\arguments{
    \item{x}{A matrix or matrix-like object containing cells in columns and features in the rows.}
    \item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
    \item{min.mean}{A numeric scalar specifying the filter to be applied on the average count for each feature prior to computing ranks.
    Disabled by setting to \code{NULL}.}
    \item{transposed}{A logical scalar specifying whether the output should be transposed.}
    \item{as.sparse}{A logical scalar indicating whether the output should be sparse.}
    \item{withDimnames}{A logical scalar specifying whether the output should contain the dimnames of \code{x}.}
}

\value{
A matrix of the same dimensions as \code{x}, where each column contains the centred and scaled ranks of the expression values for each cell.
If \code{transposed=TRUE}, this matrix is transposed so that rows correspond to cells.
}

\details{
Euclidean distances computed based on the output rank matrix are equivalent to distances computed from Spearman's rank correlation.
This can be used in clustering, nearest-neighbour searches, etc. as a robust alternative to Euclidean distances computed directly from \code{x}. 

If \code{as.sparse=TRUE}, the most common average rank is set to zero in the output.
This can be useful for highly sparse input data where zeroes have the same rank and are themselves returned as zeroes.
Obviously, this means that the ranks are not centred, so this will have to be done manually prior to any downstream distance calculations.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{quickCluster}}
}

\examples{
set.seed(100)
popsize <- 200
ngenes <- 100
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs, size=1), ncol=popsize, byrow=TRUE)

rout <- scaledColRanks(counts, transposed=TRUE)

# For use in clustering:
d <- dist(rout)
table(cutree(hclust(d), 4))

g <- buildSNNGraph(rout, transposed=TRUE)
table(igraph::cluster_walktrap(g)$membership)
}
