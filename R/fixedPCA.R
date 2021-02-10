#' PCA with a fixed number of components
#'
#' Perform a PCA where the desired number of components is known ahead of time.
#'
#' @param x A \linkS4class{SingleCellExperiment} object containing a log-expression amtrix.
#' @inheritParams denoisePCA
#' @param rank Integer scalar specifying the number of components.
#' @param preserve.shape Logical scalar indicating whether or not the output should be subsetted to \code{subset.row}.
#' Only used if \code{subset.row} is not \code{NULL}.
#'
#' @return 
#' A modified \code{x} with:
#' \itemize{
#' \item the PC results stored in the \code{\link{reducedDims}} as a \code{"PCA"} entry, if \code{type="pca"}.
#' The attributes contain the rotation matrix and 
#' \item a low-rank approximation stored as a new \code{"lowrank"} assay, if \code{type="lowrank"}.
#' This is represented as a \linkS4class{LowRankMatrix}.
#' }
#'
#' If \code{preserve.shape=TRUE}, the output always has the same number of rows as \code{x}.
#' Otherwise, the output is subsetted by any non-\code{NULL} value of \code{subset.row}.
#'
#' @details
#' In theory, there is an optimal number of components for any given application,
#' but in practice, the criterion for the optimum is difficult to define.
#' As a result, it is often satisfactory to take an \emph{a priori}-defined \dQuote{reasonable} number of PCs for downstream analyses.
#' A good rule of thumb is to set this to the upper bound on the expected number of subpopulations in the dataset
#' (see the reasoning in \code{\link{getClusteredPCs}}.
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
fixedPCA <- function(x, rank=50, value=c("pca", "lowrank"), subset.row=NULL, preserve.shape=TRUE, assay.type="logcounts", BSPARAM=bsparam(), BPPARAM=SerialParam()) {
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
