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
#' The attributes contain the rotation matrix, the variance explained and the percentage of variance explained.
#' (Note that the last may not sum to 100\% if \code{max.rank} is smaller than the total number of PCs.)
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
#' We can use \code{subset.row} to perform the PCA on a subset of genes.
#' This is typically used to subset to HVGs to reduce computational time and increase the signal-to-noise ratio of downstream analyses.
#' If \code{preserve.shape=TRUE}, the rotation matrix is extrapolated to include loadings for \dQuote{unselected} genes, i.e., not in \code{subset.row}.
#' This is done by projecting their expression profiles into the low-dimensional space defined by the SVD on the selected genes.
#' By doing so, we ensure that the output always has the same number of rows as \code{x} such that any \code{value="lowrank"} can fit into the assays.
#'
#' Otherwise, if \code{preserve.shape=FALSE}, the output is subsetted by any non-\code{NULL} value of \code{subset.row}.
#' This is equivalent to the return value after calling the function on \code{x[subset.row,]}.
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
#' @importFrom beachmat realizeFileBackedMatrix
#' @importFrom DelayedMatrixStats colVars
fixedPCA <- function(x, rank=50, value=c("pca", "lowrank"), subset.row, preserve.shape=TRUE, assay.type="logcounts", name=NULL, BSPARAM=bsparam(), BPPARAM=SerialParam()) {
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    original <- x
    x <- assay(x, assay.type)

    subset.row <- .process_subset_for_pca(subset.row, x)
    y <- t(x[subset.row,,drop=FALSE])
    y <- realizeFileBackedMatrix(y)

    svd.out <- .centered_SVD(y, rank, keep.left=TRUE, keep.right=TRUE,
        BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    var.exp <- svd.out$d^2 / (nrow(y) - 1)

    total.var <- sum(colVars(y))

    pcs <- list(
        components=.svd_to_pca(svd.out, rank), 
        rotation=.svd_to_rot(svd.out, rank, x, subset.row, fill.missing=preserve.shape),
        var.explained=var.exp,
        percent.var=var.exp/total.var*100
    )

    if (!preserve.shape) {
        original <- original[subset.row,]
    }

    value <- match.arg(value)
    .pca_to_output(original, pcs, value=value, name=name)
}
