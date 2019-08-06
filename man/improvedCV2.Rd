\name{improvedCV2}
\alias{improvedCV2}
\alias{improvedCV2,ANY-method}
\alias{improvedCV2,SingleCellExperiment-method}

\title{Stably model the technical coefficient of variation}
\description{
Model the technical coefficient of variation as a function of the mean, and determine the significance of highly variable genes.

This has been deprecated in favor of the \code{\link{modelGeneCV2}} function.
}

\usage{
\S4method{improvedCV2}{ANY}(x, is.spike, sf.cell=NULL, sf.spike=NULL, log.prior=NULL, 
    use.spikes=FALSE, top.prop=0.01, max.iter=50)

\S4method{improvedCV2}{SingleCellExperiment}(x, spike.type=NULL, ..., assay.type="logcounts", 
    logged=NULL, normalized=NULL)
}

\arguments{
\item{x}{
    A numeric matrix of counts, normalized counts or normalized log-expression values, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SingleCellExperiment object that contains such values.
}
\item{is.spike}{A vector indicating which rows of \code{x} correspond to spike-in transcripts.}
\item{sf.cell}{A numeric vector containing size factors for endogenous genes.
If this is not specified, counts for endogenous genes are assumed to already be normalized.
This is ignored if \code{log.prior!=NULL}.}
\item{sf.spike}{A numeric vector containing size factors for spike-in transcripts.
If this is not specified, counts for the spike-in transcripts are assumed to already be normalized.
This is ignored if \code{log.prior!=NULL}.}
\item{log.prior}{A numeric scalar specifying the pseudo-count added prior to log-transformation.
If this is set, \code{x} is assumed to contain normalized log-expression values, otherwise it is assumed to contain counts.}
\item{use.spikes}{A logical scalar indicating whether p-values should be returned for spike-in transcripts.}
\item{top.prop}{A numeric scalar in (0, 1), specifying the top percentage of high-abundance genes to use for estimating the initial value of the plateau parameter of the fitted trend.}
\item{max.iter}{An integer scalar specifying the maximum number of robustness iterations to use for trend fitting.}
\item{spike.type}{A character vector containing the names of the spike-in sets to use.}
\item{...}{Additional arguments to pass to \code{improvedCV2,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{logged}{A logical scalar indicating if \code{assay.type} contains log-expression values.
This is automatically determined if \code{assay.type="counts"}, \code{"logcounts"} or \code{"normcounts"}.}
\item{normalized}{A logical scalar indicating if \code{assay.type} is normalized, also automatically determined where possible.}
}

\details{
This function will estimate the squared coefficient of variation (CV2) and mean for each spike-in transcript.
Both values are log-transformed and a mean-dependent trend is fitted to the log-CV2 values. 
The trend is used to obtain the technical contribution to the CV2 for each gene.
We use a more robust approach to trend fitting that makes this function an improvement over the \code{\link{technicalCV2}} function on which it is based.

Deviations from the trend are identified by modelling the CV2 estimates for the spike-in transcripts as log-normally distributed around the fitted trend.
This accounts for sampling variance as well as any variability in the true dispersions (e.g., due to transcript-specific amplification biases).
The p-value for each gene is calculated from a one-sided Z-test on the log-CV2, using the fitted value as the mean and the robust scale estimate as the standard deviation.
A Benjamini-Hochberg adjustment is applied to correct for multiple testing.
}

\section{Trend fitting details}{
A trend of the form \eqn{y = A + B/x} is fitted to the CV2 \eqn{y} and mean expression \eqn{x}.
This represents the mean-dependent relationship that is commonly observed in RNA sequencing data.
We fit this trend using \code{\link{nls}} on the log-CV2, i.e., errors are computed as log-differences from the trend.

Several parameters are available to tune this fit:
\itemize{
\item To obtain a good initial estimate of the \eqn{A} parameter, 
we assume the CV2 has already plateaued at the top \code{top.prop} genes.
The median of these top genes is used as an initial value for \eqn{A}.
\item We reduce the effect of outliers and highly variable genes on the fit by using multiple robustness iterations.
Each gene is reweighted using tricube weights based on the residuals in log space.
Iterations are performed until convergence, for up to \code{max.iter} iterations.
}
In addition, we weight each gene in inversely proportion to the number of genes with a similar mean.
This ensures that the fit is not dominated by the large number of low-abundance genes.

For any given data set, the maximum CV2 that can be achieved is equal to the number of cells.
(This occurs when only one cell has a non-zero expression value - proof via Holder's inequality.)
Genes with CV2 values equal to the maximum are ignored during trend fitting.
This ensures that the trend is not distorted by the presence of an upper bound on CV2 values, especially at low means.
}

\section{Choice of inputs}{
For the ANY method, if \code{log.prior} is specified, \code{x} is assumed to contain log-expression values in base 2.
These are converted back to the count scale prior to calculation of the CV.
Otherwise, \code{x} is assumed to contain raw counts.
These need to be normalized with \code{sf.cell} and \code{sf.spike} prior to calculating the CV2.
Both sets of size factors are set to 1 by default if their values are not supplied. 

For \linkS4class{SingleCellExperiment} inputs, the function will try to determine the nature of the data from \code{assay.type}.
For example, \code{"counts"} are treated as raw counts and \code{"logcounts"} are treated as log2-transformed normalized values.
This can be made explicit with the \code{logged} and \code{normalized} flags.

For SingleCellExperiments with non-\code{normalized} inputs, size factors are automatically extracted via \code{\link{sizeFactors}}.
Spike-in-specific size factors for \code{spike.type} are extracted from \code{x}, if available; otherwise they are set to the size factors for the endogenous genes.
Note that the spike-in-specific factors must be the same for each set in \code{spike.type}.
}

\section{Gene selection}{
For \code{improvedCV2,ANY-method}, the rows corresponding to spike-in transcripts are specified with \code{is.spike}.
These rows will be used for trend fitting, while all other rows are treated as endogenous genes.
By default, p-values are set to \code{NA} for the spike-in transcripts, such that they do not contribute to the multiple testing correction.
This behaviour can be modified with \code{use.spikes=TRUE}, which will return p-values for all features.

For \code{improvedCV2,SingleCellExperiment-method}, transcripts from spike-in sets named in \code{spike.type} will be used for trend fitting.
If \code{spike.type=NULL}, all spike-in sets listed in \code{x} will be used.

Users can also set \code{is.spike} to \code{NA} in \code{improvedCV2,ANY-method}; or \code{spike.type} to \code{NA} in \code{improvedCV2,SingleCellExperiment-method}.
In such cases, all rows will be used for trend fitting, and (adjusted) p-values will be reported for all rows.
This should be used in cases where there are no spike-ins.
Here, the assumption is that most endogenous genes do not exhibit high biological variability and thus can be used to model decompose variation. 
}

\value{
A \linkS4class{DataFrame} is returned containing one row per row of \code{x} (including both endogenous genes and spike-in transcripts).
Each row contains the following information:
\describe{
\item{\code{mean}:}{A numeric field, containing mean (scaled) counts for all genes and transcripts.}
\item{\code{cv2}:}{A numeric field, containing CV2 values for all genes and transcripts.}
\item{\code{trend}:}{A numeric field, containing the fitted value of the trend in the CV2 values.
Note that the fitted value is reported for all genes and transcripts, but the trend is only fitted using the transcripts.}
\item{\code{ratio}:}{A numeric field, containing the ratio of \code{cv2} to \code{trend}.}
\item{\code{p.value}:}{A numeric field, containing p-values for all endogenous genes (\code{NA} for rows corresponding to spike-in transcripts).}
\item{\code{FDR}:}{A numeric field, containing adjusted p-values for all genes.}
}
}

\seealso{
\code{\link{technicalCV2}}, for the original method that this function improves on.
}

\author{
Aaron Lun
}

\examples{
# Setting up the data.
library(scRNAseq)
sce <- ZeiselBrainData()

library(scater)
sizeFactors(sce) <- librarySizeFactors(sce)
sce <- computeSpikeFactors(sce)

# Computing based on the spike-ins.
out <- improvedCV2(sce, assay.type="counts")
head(out)

plot(out$mean, out$cv2, log="xy")
points(out$mean, out$trend, col="red", pch=16, cex=0.5)

# Again, but without the spike.ins.
out2 <- improvedCV2(sce, spike.type=NA, assay.type="counts")
plot(out2$mean, out2$cv2, log="xy")
points(out2$mean, out2$trend, col="red", pch=16, cex=0.5)
}

\keyword{variance}
