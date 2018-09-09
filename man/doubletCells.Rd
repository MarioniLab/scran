\name{doubletCells}
\alias{doubletCells}
\alias{doubletCells,ANY-method}
\alias{doubletCells,SingleCellExperiment-method}

\title{Detect doublet cells}
\description{Identify potential doublet cells based on simulations of putative doublet expression profiles.}

\usage{
\S4method{doubletCells}{ANY}(x, size.factors.norm=NULL, size.factors.content=NULL,
    k=50, subset.row=NULL, niters=max(10000, ncol(x)), block=10000, 
    d=50, approximate=FALSE, irlba.args=list(), force.match=FALSE, 
    force.k=20, force.ndist=3, BPPARAM=SerialParam())

\S4method{doubletCells}{SingleCellExperiment}(x, size.factors.norm=NA, ..., subset.row=NULL, 
    assay.type="counts", get.spikes=FALSE)
}

\arguments{
\item{x}{
A numeric matrix-like object of count values, where each column corresponds to a cell and each row corresponds to an endogenous gene.
Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{size.factors.norm}{A numeric vector of size factors for normalization of \code{x} prior to PCA and distance calculations.
If \code{NULL}, defaults to the column sums of \code{x}.

For the SingleCellExperiment method, this may be \code{NA}, in which case \code{sizeFactors(x)} is used instead.}
\item{size.factors.content}{A numeric vector of size factors for RNA content normalization of \code{x} prior to simulating doublets.}
\item{k}{An integer scalar specifying the number of nearest neighbours to use to determine the bandwidth for density calculations.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{niters}{An integer scalar specifying how many simulated doublets should be generated.}
\item{block}{An integer scalar controlling the rate of doublet generation, to keep memory usage low.}
\item{d}{An integer scalar specifying the number of components to retain after the PCA.}
\item{approximate}{A logical scalar indicating whether \code{\link[irlba]{irlba}} should be used to perform the initial PCA.}
\item{irlba.args}{A list of arguments to pass to \code{\link[irlba]{irlba}} when \code{approximate=TRUE}.}
\item{force.match}{A logical scalar indicating whether remapping of simulated doublets to original cells should be performed.}
\item{force.k}{An integer scalar specifying the number of neighbours to use for remapping if \code{force.match=TRUE}.}
\item{force.ndist}{A numeric scalar specifying the bandwidth for remapping if \code{force.match=TRUE}.}
\item{BPPARAM}{A BiocParallelParam object specifying whether the nearest-neighbour searches should be parallelized.}
\item{...}{Additional arguments to pass to the ANY method.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}
\item{get.spikes}{A logical scalar specifying whether decomposition should be performed for spike-ins.}
}

\value{
A numeric vector of doublet scores for each cell in \code{x}.
}

\details{
This function simulates doublets by adding the count vectors for two randomly chosen cells in \code{x}.
For each cell, we compute the density of simulated doublets and compare it to the density of original cells.
Genuine doublets should have a high density of simulated doublets relative to the density of its neighbourhood.
Thus, the doublet score for each cell is defined as the ratio of densities of simulated doublets to the (squared) density of the original cells.

Densities are calculated in low-dimensional space after a PCA on the log-normalized expression matrix of \code{x}.
Simulated doublets are projected into the low-dimensional space using the rotation vectors computed from the original cells.
A tricube kernel is used to compute the density around each cell.
The bandwidth of the kernel is set to the median distance to the \code{k} nearest neighbour across all cells.

The two size factor arguments have different roles:
\itemize{
\item \code{size.factors.norm} contains the size factors to be used for normalization prior to PCA and distance calculations.
This can be set to ensure that the low-dimensional space is consistent with that in the rest of the analysis.
\item \code{size.factors.content} is much more important, and represents the size factors that preserve RNA content differences.
This is usually computed from spike-in RNA and ensures that the simulated doublets have the correct ratio of contributions from the original cells.
}
It is possible to set both of these arguments, as they will not interfere with each other.

If \code{force.match=TRUE}, simulated doublets will be remapped to the nearest neighbours in the original data.
This is done by taking the (tricube-weighted) average of the PC scores for the \code{force.k} nearest neighbors.
The tricube bandwidth for remapping is chosen by taking the median distance and multiplying it by \code{force.ndist}, to protect against later neighbours that might be outliers.
The aim is to adjust for unknown differences in RNA content that would cause the simulated doublets to be systematically displaced from their true locations.
However, it may also result in spuriously high scores for single cells that happen to be close to a cluster of simulated doublets.
}

\author{
Aaron Lun
}

\examples{
# Mocking up an example.
ngenes <- 100
mu1 <- 2^rexp(ngenes)
mu2 <- 2^rnorm(ngenes)

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes)
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes)
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes)

counts <- cbind(counts.1, counts.2, counts.m)
clusters <- rep(1:3, c(ncol(counts.1), ncol(counts.2), ncol(counts.m)))

# Find potential doublets...
scores <- doubletCells(counts)
boxplot(split(scores, clusters))
}
