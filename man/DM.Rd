\name{Distance-to-median}
\alias{DM}

\title{Compute the distance-to-median statistic}
\description{Compute the distance-to-median statistic for the CV2 residuals of all genes}

\usage{
DM(mean, cv2, win.size=50)
}

\arguments{
\item{mean}{A numeric vector of average counts for each gene.}
\item{cv2}{A numeric vector of squared coefficients of variation for each gene.}
\item{win.size}{An integer scalar specifying the window size for median-based smoothing.}
}

\details{
This function will compute the distance-to-median (DM) statistic described by Kolodziejczyk et al. (2015).
Briefly, a median-based trend is fitted to the log-transformed \code{cv2} against the log-transformed \code{mean}.
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
ngenes <- 10000
means <- exp(runif(ngenes, 2, 8))
cv2 <- rgamma(ngenes, 5, 5)
dm.stat <- DM(means, cv2)
}

\references{
Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015). 
Single cell RNA-sequencing of pluripotent states unlocks modular transcriptional variation. 
\emph{Cell Stem Cell} 17(4), 471--85.
}

\keyword{variance}
