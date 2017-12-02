\name{Distance-to-median}
\alias{DM}

\title{Compute the distance-to-median statistic}
\description{Compute the distance-to-median statistic for the CV2 residuals of all genes}

\usage{
DM(mean, cv2, win.size=51)
}

\arguments{
\item{mean}{A numeric vector of average counts for each gene.}
\item{cv2}{A numeric vector of squared coefficients of variation for each gene.}
\item{win.size}{An integer scalar specifying the window size for median-based smoothing.
This should be odd or will be incremented by 1.} 
}

\details{
This function will compute the distance-to-median (DM) statistic described by Kolodziejczyk et al. (2015).
Briefly, a median-based trend is fitted to the log-transformed \code{cv2} against the log-transformed \code{mean} using \code{\link{runmed}}.
The DM is defined as the residual from the trend for each gene.
This statistic is a measure of the relative variability of each gene, after accounting for the empirical mean-variance relationship.
Highly variable genes can then be identified as those with high DM values.
}

\value{
A numeric vector of DM statistics for all genes.
}

\author{
Jong Kyoung Kim,
with modifications by Aaron Lun
}

\examples{
# Mocking up some data
ngenes <- 1000
ncells <- 100
gene.means <- 2^runif(ngenes, 0, 10)
dispersions <- 1/gene.means + 0.2
counts <- matrix(rnbinom(ngenes*ncells, mu=gene.means, size=1/dispersions), nrow=ngenes)

# Computing the DM.
means <- rowMeans(counts)
cv2 <- apply(counts, 1, var)/means^2
dm.stat <- DM(means, cv2)
head(dm.stat)
}

\references{
Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015). 
Single cell RNA-sequencing of pluripotent states unlocks modular transcriptional variation. 
\emph{Cell Stem Cell} 17(4), 471--85.
}

\keyword{variance}
