\name{trendVar}
\alias{trendVar}
\alias{trendVar,ANY-method}
\alias{trendVar,SingleCellExperiment-method}

\title{Fit a variance trend}
\description{Fit a mean-dependent trend to the gene-specific variances in single-cell RNA-seq data.}

\usage{
\S4method{trendVar}{ANY}(x, method=c("loess", "spline"), parametric=FALSE, 
    loess.args=list(), spline.args=list(), nls.args=list(), block=NULL, 
    design=NULL, weighted=TRUE, min.mean=0.1, subset.row=NULL, 
    BPPARAM=SerialParam()) 

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
\item{spline.args}{A named list of arguments to pass to \code{\link[aroma.light]{robustSmoothSpline}} when \code{method="spline"}.}
\item{nls.args}{A named list of arguments to pass to \code{\link{nls}} when \code{parametric=TRUE}.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{design}{A numeric matrix describing the uninteresting factors contributing to expression in each cell.
Alternatively, a single factor can be supplied for one-way layouts.}
\item{weighted}{A logical scalar indicated whether weighted trend fitting should be performed when \code{block!=NULL}.}
\item{min.mean}{A numeric scalar specifying the minimum mean log-expression in order for a gene to be used for trend fitting.}
\item{subset.row}{See \code{?"\link{scran-gene-selection}"}.}
\item{BPPARAM}{A BiocParallelParam object indicating whether and how parallelization should be performed across genes.}
\item{...}{Additional arguments to pass to \code{trendVar,ANY-method}.}
\item{assay.type}{A string specifying which assay values in \code{x} to use.}
\item{use.spikes}{A logical scalar specifying whether the trend should be fitted to variances for spike-in transcripts or endogenous genes.}
}

\details{
This function fits an abundance-dependent trend to the variance of the log-normalized expression for the spike-in transcripts.
For SingleCellExperiment objects, these expression values are computed by \code{\link[=normalize,SingleCellExperiment-method]{normalize}} after setting the size factors, 
e.g., with \code{\link{computeSpikeFactors}}.
Log-transformed values are used as these are more robust to genes/transcripts with strong expression in only one or two outlier cells.
It also allows the fitted trend to be applied in downstream procedures that use log-transformed counts.

The mean and variance of the normalized log-counts is calculated for each spike-in transcript, and a trend is fitted to the variance against the mean for all transcripts.
The fitted value of this trend represents technical variability due to sequencing, drop-outs during capture, etc. at a given mean.
This assumes that a constant amount of spike-in RNA was added to each cell, such that any differences in observed expression are purely due to measurement error.
Variance decomposition to biological and technical components for endogenous genes can then be performed later with \code{\link{decomposeVar}}.
}

\section{Trend fitting options}{
If \code{parametric=FALSE}, smoothing is performed directly on the log-variances.
This is the default as it provides the most stable performance on arbitrary mean-variance relationships.

If \code{parametric=TRUE}, a non-linear curve of the form
\deqn{y = \frac{ax}{x^n + b}}{y = ax/(x^n + b)}
is fitted to the variances against the means using \code{\link{nls}}.
Starting values and the number of iterations are automatically set if not explicitly specified in \code{nls.args}.
A smoothing algorithm is then applied to the log-ratios of the variance to the fitted value for each gene.
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
Spline fitting is also faster for very large numbers of points, which may occur in large data sets with many genes and/or many levels of \code{block}.

The \code{trendVar} function will produce an output \code{trend} function with which fitted values can be computed.
When extrapolating to values below the smallest observed mean, the output function will approach zero.
When extrapolating to values above the largest observed mean, the output function will be set to the fitted value of the trend at the largest mean.
}

\section{Handling uninteresting factors of variation}{
There are three approaches to handling unwanted factors of variation.
The simplest approach is to use a design matrix containing the uninteresting factors can be specified in \code{design}.
This will fit a linear model to the log-expression values for each gene, yielding an estimate for the residual variance.
The trend is then fitted to the residual variance against the mean for each spike-in transcripts.

Another approach is to use \code{block}, where all cells in each level of the blocking factor are treated as a separate group.
Means and variances are estimated within each group and the resulting sets of means/variances are pooled across all groups.
The trend is then fitted to the pooled observations, where observations from different levels are weighted according to the residual d.f. used for variance estimation.
This effectively multiplies the number of points by the number of levels in \code{block}.
If both \code{block} and \code{design} are specified, \code{block} will take priority and \code{design} will be ignored.

The final approach is to subset the data set for each level of the blocking factor and model the variances within each block separately.
This is done using the \code{\link{multiBlockVar}} function, see its documentation for more details.
Separate modelling and trend fitting is the most correct approach if there are systematic differences in the size factors (spike-in or endogenous) between levels.
With the other two methods, such differences would be normalized out in the full log-expression matrix, preventing proper estimation of the level-specific abundance.

Assuming there are no differences in the size factors between levels, we suggest using \code{block} wherever possible instead of \code{design}.
This is because the use of \code{block} preserves differences in the means/variances between levels of the factor. 
In contrast, using \code{design} will effectively compute an average mean/variance.
This may yield an inaccurate representation of the trend, as the fitted value at an average mean may not be equal to the average variance for non-linear trends.
Nonetheless, we still support \code{design} as it can accommodate additive models, whereas \code{block} only handles one-way layouts.

% This function is separated from decomposeVar in order to allow for different design matrices/cells to be used for trend fitting and variance estimation.
% For example, if you have multiple groups, you can fit the trend to the spike-ins for all groups, but estimate the variance individually for each group.
}

\section{Additional notes on row selection}{
The selection of spike-in transcripts can be adjusted in \code{trendVar,SingleCellExperiment-method} using the \code{use.spikes} method. 
\itemize{
\item By default, \code{use.spikes=TRUE} which means that only rows labelled as spike-ins with \code{isSpike(x)} will be used.
An error will be raised if no rows are labelled as spike-in transcripts.
\item If \code{use.spikes=FALSE}, only the rows \emph{not} labelled as spike-in transcripts will be used.
\item If \code{use.spikes=NA}, every row will be used for trend fitting, regardless of whether it corresponds to a spike-in transcript or not.
}

If \code{use.spikes=FALSE}, this implies that \code{trendVar} will be applied to the endogenous genes in the SingleCellExperiment object.
For \code{trendVar,ANY-method}, it is equivalent to manually supplying a matrix of normalized expression for endogenous genes.
This assumes that most genes exhibit technical variation and little biological variation, e.g., in a homogeneous population.

Low-abundance genes with mean log-expression below \code{min.mean} are not used in trend fitting, to preserve the sensitivity of span-based smoothers at moderate-to-high abundances.
It also protects against discreteness, which can interfere with estimation of the variability of the variance estimates and accurate scaling of the trend.
The default threshold is chosen based on the point at which discreteness is observed in variance estimates from Poisson-distributed counts.
For heterogeneous droplet data, a lower threshold of 0.001-0.01 should be used.

Users can also directly specify which rows to use with \code{subset.row}.
All of these parameters - including \code{min.mean} and \code{use.spikes} - interact with each other sensibly, by taking the intersection of rows that fulfill all criteria.
For example, if \code{subset.row} is specified and \code{use.spikes=TRUE}, rows are only used if they are both spike-in transcripts \emph{and} in \code{subset.row}.
Otherwise, if \code{use.spikes=FALSE}, only rows in \code{subset.row} that are \emph{not} spike-in transcripts are used.

Note that features are not used under any circumstance if they have a variance of zero.
These are not compatible with trend fitting on the log-scale.
In practice, it should be effectively impossible to obtain variances of zero for any expressed gene, though some care may be required with simulated data sets.
Any \code{NA} variances are similarly ignored.

Note that, depending on the output of the various row selection parameters, the number of used rows may be much less than the number of rows of \code{x}.
In pathological cases, this may result in an error stating that insufficient points are available for the trend points.
}

\section{Warning on size factor centring}{ 
If \code{assay.type="logcounts"}, \code{trendVar,SingleCellExperiment-method} will attempt to determine if the expression values were computed from counts via \code{\link[scater]{normalize}}.
If so, a warning will be issued if the size factors are not centred at unity.
This is because different size factors are typically used for endogenous genes and spike-in transcripts.
If these size factor sets are not centred at the same value, there will be systematic differences in abundance between these features.
This precludes the use of a spike-in fitted trend with abundances for endogenous genes in \code{\link{decomposeVar}}.

For other expression values and in \code{trendVar,ANY-method}, the onus is on the user to ensure that normalization 
(i) does not introduce differences in abundance between spike-in and endogenous features, while
(ii) preserving differences in abundance within the set of endogenous or spike-in features.
In short, the scaling factors used to normalize each feature should have the same mean across all cells.
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
\item{\code{mean}:}{A numeric vector of mean log-expression values for all spike-in transcripts, if \code{block=NULL}.
Otherwise, a numeric matrix of means where each row corresponds to a spike-in and each column corresponds to a level of \code{block}.}
\item{\code{var}:}{A numeric vector of the variances of log-expression values for all spike-in transcripts, if \code{block=NULL}.
Otherwise, a numeric matrix of variances where each row corresponds to a spike-in and each column corresponds to a level of \code{block}.}
\item{\code{resid.df}:}{An integer scalar specifying the residual d.f. used for variance estimation of each spike-in transcript, if \code{block=NULL}.
Otherwise, a integer vector where each entry specifies the residual d.f. used in each level of \code{block}.}
\item{\code{block}:}{A factor identical to the input \code{block}, only returned if it was not \code{NULL}.}
\item{\code{design}:}{A numeric matrix (or factor) identical to the input \code{design}, only returned if it was not \code{NULL} and \code{block=NULL}.}
\item{\code{trend}:}{A function that returns the fitted value of the trend at any mean.}
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

