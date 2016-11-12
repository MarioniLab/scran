\name{sandbag}
\alias{sandbag}
\alias{sandbag,matrix-method}
\alias{sandbag,SCESet-method}

\title{Cell cycle phase training}
\description{Use gene expression data to train a classifier for cell cycle phase.}

\usage{
\S4method{sandbag}{matrix}(x, is.G1, is.S, is.G2M, gene.names=rownames(x), 
    fraction=0.5, subset.row=NULL)

\S4method{sandbag}{SCESet}(x, is.G1, is.S, is.G2M, subset.row=NULL, ..., 
    assay="counts", get.spikes=FALSE)
}

\arguments{
\item{x}{
    A numeric matrix of gene expression values where rows are genes and columns are cells.
    Alternatively, a SCESet object containing such a matrix.
}
\item{is.G1, is.S, is.G2M}{A vector indicating which cells are in each phase of the cell cycle.}
\item{gene.names}{A character vector of gene names.}
\item{fraction}{A numeric scalar specifying the minimum fraction to define a marker gene pair.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to \code{sandbag,matrix-method}.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{get.spikes}{A logical specifying whether spike-in transcripts should be used.}
}

\details{
This function implements the training step of the pair-based prediction method described by Scialdone et al. (2015).
Pairs of genes (A, B) are identified from a training data set where in each pair,
    the fraction of cells in phase G1 with expression of A > B (based on expression values in \code{training.data}) 
    and the fraction with B > A in each other phase exceeds \code{fraction}.
These pairs are defined as the marker pairs for G1.
This is repeated for each phase to obtain a separate marker pair set.

Pre-defined sets of marker pairs are provided for mouse and human (see Examples).
The mouse set was generated as described by Scialdone et al. (2015), while the human training set was generated with data from Leng et al. (2015).
Classification from test data can be performed using the \code{\link{cyclone}} function.
For each cell, this involves comparing expression values between genes in each marker pair. 
The cell is then assigned to the phase that is consistent with the direction of the difference in expression in the majority of pairs.

For \code{sandbag,SCESet-method}, the matrix of counts is used but can be replaced with expression values by setting \code{assays}.
By default, \code{get.spikes=FALSE} which means that any rows corresponding to spike-in transcripts will not be considered when picking markers.
This is because the amount of spike-in RNA added will vary between experiments and will not be a robust predictor.
Nonetheless, if all rows are required, users can set \code{get.spikes=TRUE}.
Users can also manually select which rows to use via \code{subset.row}, which will override any setting of \code{get.spikes}.
}

\value{
A named list of data.frames, where each data frame corresponds to a cell cycle phase and contains the names of the genes in each marker pair.
}

\author{
Antonio Scialdone,
with modifications by Aaron Lun
}

\seealso{
\code{\link{cyclone}}
}

\examples{
ncells <- 50
ngenes <- 20
training <- matrix(rnorm(ncells*ngenes), ncol=ncells)
rownames(training) <- paste0("X", seq_len(ngenes))

is.G1 <- 1:20
is.S <- 21:30
is.G2M <- 31:50
out <- sandbag(training, is.G1, is.S, is.G2M) 
str(out)

# Getting pre-trained marker sets
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
}

\references{
Scialdone A, Natarajana KN, Saraiva LR et al. (2015). 
Computational assignment of cell-cycle stage from single-cell transcriptome data.
\emph{Methods} 85:54--61

Leng N, Chu LF, Barry C et al. (2015).
Oscope identifies oscillatory genes in unsynchronized single-cell RNA-seq experiments.
\emph{Nat. Methods} 12:947--50
}

\keyword{clustering}
