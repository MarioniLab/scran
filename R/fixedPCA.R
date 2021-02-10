#' PCA with a fixed number of components
#'
#' Perform a PCA where the desired number of components is known ahead of time.
#'
#' @param x A \linkS4class{SingleCellExperiment} object containing a log-expression amtrix.
#' @inheritParams denoisePCA
#' @param rank Integer scalar specifying the number of components.
#'
#' @return 
#' A modified \code{x} with:
#' \itemize{
#' \item the PC results stored in the \code{\link{reducedDims}} as a \code{"PCA"} entry, if \code{type="pca"}.
#' The attributes contain the rotation matrix and the percentage of variance explained.
#' \item a low-rank approximation stored as a new \code{"lowrank"} assay, if \code{type="lowrank"}.
#' This is represented as a \linkS4class{LowRankMatrix}.
#' }
#'
#' @details
#' In theory, there is an optimal number of components for any given application,
#' but in practice, the criterion for the optimum is difficult to define.
#' As a result, it is often satisfactory to take an \emph{a priori}-defined \dQuote{reasonable} number of PCs for downstream analyses.
#' A good rule of thumb is to set this to the upper bound on the expected number of subpopulations in the dataset
#' (see the reasoning in \code{\link{getClusteredPCs}}.
#' 
#' If \code{preserve.shape=TRUE}, the rotation matrix is extrapolated to include loadings for \dQuote{unselected} genes, i.e., not in \code{subset.row}.
#' This is done by projecting their expression profiles into the low-dimensional space defined by the SVD on the selected genes.
#' By doing so, we ensure that the output always has the same number of rows as \code{x} such that any \code{value="lowrank"} can fit into the assays.
#' Otherwise, the output is subsetted by any non-\code{NULL} value of \code{subset.row}.
#'
#' @author Aaron Lun
#'
#' @seealso 
#' \code{\link{denoisePCA}}, where the number of PCs is automatically chosen.
#'
#' \code{\link{getClusteredPCs}}, another method to choose the number of PCs.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' # Modelling the variance:
#' var.stats <- modelGeneVar(sce)
#' hvgs <- getTopHVGs(var.stats, n=1000)
#'
#' # Defaults to pulling out the top 50 PCs.
#' set.seed(1000)
#' sce <- fixedPCA(sce, subset.row=hvgs)
#' reducedDimNames(sce)
#'
#' # Get the percentage of variance explained. 
#' attr(reducedDim(sce), "percentVar")
#' 
#' @export
#' @importFrom BiocSingular bsparam
#' @importFrom BiocParallel SerialParam bpstop bpstart
#' @importFrom SummarizedExperiment assay
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom Matrix t 
fixedPCA <- function(x, rank=50, value=c("pca", "lowrank"), subset.row=NULL, preserve.shape=(value=="lowrank"), assay.type="logcounts", BSPARAM=bsparam(), BPPARAM=SerialParam()) {
    if (!.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    original <- x
    x <- assay(x, assay.type)

    subset.row <- .subset2index(subset.row, x, byrow=TRUE)
    stats <- .compute_mean_var(x, BPPARAM=BPPARAM, subset.row=subset.row, 
        design=NULL, block.FUN=compute_blocked_stats_none, block=NULL)
    total.var <- sum(stats$vars)

    y <- x[subset.row,,drop=FALSE]
    svd.out <- .centered_SVD(t(y), rank, keep.left=TRUE, keep.right=TRUE,
        BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    var.exp <- svd.out$d^2 / (ncol(y) - 1)

    value <- match.arg(value)
    force(preserve.shape) # must happen after 'value', for the defaults to work properly.

    pcs <- list(
        components=.svd_to_pca(svd.out, rank), 
        rotation=.svd_to_rot(svd.out, rank, x, subset.row, fill.missing=preserve.shape),
        percent.var=var.exp/total.var*100
    )

    if (!preserve.shape) {
        original <- original[subset.row,]
    }
    .pca_to_output(original, pcs, value=value)
}
