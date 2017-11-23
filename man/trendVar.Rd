\name{trendVar}
\alias{trendVar}
\alias{trendVar,ANY-method}
\alias{trendVar,SingleCellExperiment-method}

\title{Fit a variance trend}
\description{Fit a mean-dependent trend to the gene-specific variances in single-cell RNA-seq data.}

\usage{
\S4method{trendVar}{ANY}(x, method=c("loess", "spline"), parametric=FALSE, 
    loess.args=list(), spline.args=list(), nls.args=list(),
    span=NULL, family=NULL, degree=NULL, df=NULL, start=NULL, 
    min.mean=0.1, design=NULL, subset.row=NULL) 

\S4method{trendVar}{SingleCellExperiment}(x, subset.row=NULL, ..., assay.type="logcounts", use.spikes=TRUE)
}

\arguments{
\item{x}{
    A numeric matrix-like object of normalized log-expression values, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SingleCellExperiment object that contains such values.
}
\item{method}{A string specifying the algorithm to use for smooth trend fitting.}
\item{parametric}{A logical scalar indicating whether a parametric curve should be fitted prior to smoothing.}
\item{loess.args}{A named list of arguments to pass to \code{\link{loess}} when \code{method="loess"}.}
\item{spline.args}{A named list of arguments to pass to \code{\link[stats]{smooth.spline}} when \code{method="spline"}.}
\item{nls.args}{A named list of arguments to pass to \code{\link{nls}} when \code{parametric=TRUE}.}
\item{span, family, degree}{Deprecated arguments to pass to \code{\link{loess}}.}
\item{df}{Deprecated argument to pass to \code{\link{ns}}.}
\item{start}{Deprecated argument to pass to \code{\link{nls}}.}
\item{min.mean}{A numeric scalar specifying the minimum mean log-expression in order for a gene to be used for trend fitting.}
\item{design}{A numeric matrix describing the uninteresting factors contributing to expression in each cell.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to \code{trendVar,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{use.spikes}{A logical scalar specifying whether the trend should be fitted to variances for spike-in transcripts or endogenous genes.}
}

\details{
This function fits an abundance-dependent trend to the variance of the log-normalized expression for the spike-in transcripts.
For SingleCellExperiment objects, these expression values can be computed by \code{\link{normalize}} after setting the size factors, e.g., with \code{\link{computeSpikeFactors}}.
Log-transformed values are used as these are more robust to genes/transcripts with strong expression in only one or two outlier cells.

The mean and variance of the normalized log-counts is calculated for each spike-in transcript, and a trend is fitted to the variance against the mean for all transcripts.
The fitted value of this trend represents technical variability due to sequencing, drop-outs during capture, etc.
Variance decomposition to biological and technical components for endogenous genes can then be performed later with \code{\link{decomposeVar}}.

The design matrix can be set if there are factors that should be blocked, e.g., batch effects, known (and uninteresting) clusters.
Otherwise, it will default to an all-ones matrix, effectively treating all cells as part of the same group.

% This function is separated from decomposeVar in order to allow for different design matrices/cells to be used for trend fitting and variance estimation.
% For example, if you have multiple groups, you can fit the trend to the spike-ins for all groups, but estimate the variance individually for each group.
}

\section{Trend fitting options}{
If \code{parametric=FALSE}, smoothing is performed directly on the log-variances.
This is the default as it provides the most stable performance on arbitrary mean-variance relationships.

If \code{parametric=TRUE}, a non-linear curve of the form
\deqn{y = \frac{ax}{x^n + b}}{y = ax/(x^n + b)}
is fitted to the variances against the means using \code{\link{nls}}.
Starting values and the number of iterations are automatically set if not explicitly specified in \code{nls.args}.
A smoothing algorihtm is then applied to the log-ratios of the variance to the fitted value for each gene.
The aim is to use the parametric curve to reduce the sharpness of the expected mean-variance relationship[for easier smoothing.
Conversely, the parametric form is not exact, so the smoothers will model any remaining trends in the residuals.

The \code{method} argument specifies the smoothing algorithm to be applied on the log-ratios/variances.
By default, a robust loess curve is used for trend fitting via \code{\link{loess}}.
This provides a fairly flexible fit while protecting against genes with very large or very small variances.
Arguments to \code{\link{loess}} are specified with \code{loess.args}, with defaults of \code{span=0.3}, \code{family="symmetric"} and \code{degree=1} unless otherwise specified. 
Some experimentation with these parameters may be required to obtain satisfactory results.

If \code{method="spline"}, smoothing will instead be performed using the \code{\link{smooth.spline}} function 
Arguments are specified with \code{spline.args}, with a default degrees of freedom of \code{df=4} unless otherwise specified.
Splines can be more effective than loess at capturing smooth curves with strong non-linear gradients.

The \code{trendVar} function will produce an output \code{trend} function with which fitted values can be computed.
When extrapolating to values below the smallest observed mean, the output function will approach zero.
When extrapolating to values above the largest observed mean, the output function will be set to the fitted value of the trend at the largest mean.
}

\section{Additional notes on row selection}{
Spike-in transcripts can be selected in \code{trendVar,SingleCellExperiment-method} using the \code{use.spikes} method. 
\itemize{
\item By default, \code{use.spikes=TRUE} which means that only rows labelled as spike-ins with \code{isSpike(x)} will be used.
\item If \code{use.spikes=FALSE}, only the rows \emph{not} labelled as spike-ins will be used.
\item If \code{use.spikes=NA}, every row will be used for trend fitting, regardless of whether it corresponds to a spike-in transcript or not.
}

If \code{use.spikes=FALSE}, this implies that \code{trendVar} will be applied to the endogenous genes in the SingleCellExperiment object.
For \code{trendVar,ANY-method}, it is equivalent to manually supplying a matrix of normalized expression for endogenous genes.
This assumes that most genes exhibit technical variation and little biological variation, e.g., in a homogeneous population.

Users can also directly specify which rows to use with \code{subset.row}, which will interact as expected with any non-\code{NA} value of \code{use.spikes}.
If \code{subset.row} is specified and \code{use.spikes=TRUE}, only the spike-in transcripts in \code{subset.row} are used.
Otherwise, if \code{use.spikes=FALSE}, only the non-spike in transcripts in \code{subset.row} are used.

Low-abundance genes with mean log-expression below \code{min.mean} are not used in trend fitting, to preserve the sensitivity of span-based smoothers at moderate-to-high abundances.
It also protects against discreteness, which can interfere with estimation of the variability of the variance estimates and accurate scaling of the trend.
The default threshold is chosen based on the point at which discreteness is observed in variance estimates from Poisson-distributed counts.
For heterogeneous droplet data, a lower threshold of 0.001-0.01 should be used.

Obviously, the usefulness of the trend is dependent on the quality of the features to which it is fitted.
For example, the trend will not provide accurate estimates of the technical component if the coverage of all spike-ins is lower than that of endogenous genes.
}

\section{Warning on size factor centring}{ 
If \code{assay.type="logcounts"}, \code{trendVar,SingleCellExperiment-method} will attempt to determine if the expression values were computed from counts via \code{\link[scater]{normalize}}.
If so, a warning will be issued if the size factors are not centred at unity.
This is because different size factors are typically used for endogenous genes and spike-in transcripts.
If these size factor sets are not centred at the same value, there will be systematic differences in abundance between these features.
This precludes the use of a spike-in fitted trend with abundances for endogenous genes in \code{\link{decomposeVar}}.

For other expression values and in \code{trendVar,ANY-method}, the onus is on the user to ensure that normalization preserves differences in abundance.
In other words, the scaling factors used to normalize each feature should have the same mean.
This ensures that spurious differences in abundance are not introduced by the normalization process.

% Forcibly recomputing the "logcounts" values within trendVar() is arguably safer, as normalize() enforces centring.
% However, it's slightly expensive, hidden from the user and restricted to this function.
% Better to request a re-calculation on the user's side, such that the recalculated values are consistently used here and elsewhere.
%
% The attempted check should only fail in pathological cases:
% - false positives when counts are added to "x" after computing "logcounts" (fairly unusual to not put them in as raw data first)
% - false positives when normalize() is run with centre_size_factors=FALSE (don't know why you would do that)
% - false negatives when counts are removed from "x" after computing "logcounts" (this is pretty stupid)
% - false negatives when size factors are modified by centring but "logcounts" is not re-computed (arguably a problem for the entire workflow)
%
% The other options would involve enforcing validity checks on "logcounts" within the SingleCellExperiment object itself.
% However, this is probably too obtrusive for little gain; only this function requires consistency in the absolute values of the size factors.
}

\value{
A named list is returned, containing:
\describe{
\item{\code{mean}:}{A numeric vector of mean log-CPMs for all spike-in transcripts.}
\item{\code{var}:}{A numeric vector of the variances of log-CPMs for all spike-in transcripts.}
\item{\code{trend}:}{A function that returns the fitted value of the trend at any mean log-CPM.}
\item{\code{design}:}{A numeric matrix, containing the design matrix that was used.}
\item{\code{df2}:}{A numeric scalar, specifying the second degrees of freedom for a scaled F-distribution describing the variability of variance estimates around the trend.}
}
}

\seealso{
\code{\link{nls}},
\code{\link{loess}},
\code{\link{decomposeVar}},
\code{\link{computeSpikeFactors}},
\code{\link{computeSumFactors}},
\code{\link{normalize}}
}

\author{
Aaron Lun
}

\examples{
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.

# Normalizing (gene-based factors for genes, spike-in factors for spike-ins)
y <- computeSumFactors(y) 
y <- computeSpikeFactors(y, general.use=FALSE)
y <- normalize(y)

# Fitting a trend to the spike-ins.
fit <- trendVar(y)
plot(fit$mean, fit$var)
curve(fit$trend(x), col="red", lwd=2, add=TRUE)

# Fitting a trend to the endogenous genes. 
fit.g <- trendVar(y, use.spikes=FALSE)
plot(fit.g$mean, fit.g$var)
curve(fit.g$trend(x), col="red", lwd=2, add=TRUE)
}

\references{
Lun ATL, McCarthy DJ and Marioni JC (2016). 
A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.
\emph{F1000Res.} 5:2122

}

\keyword{variance}

