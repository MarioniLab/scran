#' Compute scaled column ranks
#'
#' Compute scaled column ranks from each cell's expression profile for distance calculations based on rank correlations.
#' 
#' @param x A numeric matrix-like object containing cells in columns and features in the rows.
#' @param subset.row A logical, integer or character scalar indicating the rows of \code{x} to use, see \code{?"\link{scran-gene-selection}"}.
#' @param min.mean A numeric scalar specifying the filter to be applied on the average normalized count for each feature prior to computing ranks.
#' Disabled by setting to \code{NULL}.
#' @param transposed A logical scalar specifying whether the output should be transposed.
#' @param as.sparse A logical scalar indicating whether the output should be sparse.
#' @param withDimnames A logical scalar specifying whether the output should contain the dimnames of \code{x}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether and how parallelization should be performed.
#' Currently only used for filtering if \code{min.mean} is not provided.
#' 
#' @return
#' A matrix of the same dimensions as \code{x}, where each column contains the centred and scaled ranks of the expression values for each cell.
#' If \code{transposed=TRUE}, this matrix is transposed so that rows correspond to cells.
#' If \code{as.sparse}, the columns are not centered to preserve sparsity.
#' 
#' @details
#' Euclidean distances computed based on the output rank matrix are equivalent to distances computed from Spearman's rank correlation.
#' This can be used in clustering, nearest-neighbour searches, etc. as a robust alternative to Euclidean distances computed directly from \code{x}. 
#' 
#' If \code{as.sparse=TRUE}, the most common average rank is set to zero in the output.
#' This can be useful for highly sparse input data where zeroes have the same rank and are themselves returned as zeroes.
#' Obviously, this means that the ranks are not centred, so this will have to be done manually prior to any downstream distance calculations.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{quickCluster}}, where this function is used.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' rout <- scaledColRanks(counts(sce), transposed=TRUE)
#' 
#' # For use in clustering:
#' d <- dist(rout)
#' table(cutree(hclust(d), 4))
#' 
#' g <- buildSNNGraph(rout, transposed=TRUE)
#' table(igraph::cluster_walktrap(g)$membership)
#' 
#' @export
#' @importFrom scuttle calculateAverage .subset2index .bpNotSharedOrUp
#' @importFrom BiocParallel SerialParam bpstart bpstop
scaledColRanks <- function(x, subset.row=NULL, min.mean=NULL, transposed=FALSE, as.sparse=FALSE, 
    withDimnames=TRUE, BPPARAM=SerialParam())
{
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    subset.row <- .subset2index(subset.row, x, byrow=TRUE)
    if (!is.null(min.mean) && all(dim(x)>0L)) {
        further.subset <- calculateAverage(x, subset_row=subset.row, BPPARAM=BPPARAM) >= min.mean
        subset.row <- subset.row[further.subset]
    }

    out <- colBlockApply(x[subset.row,,drop=FALSE], FUN=.get_scaled_ranks, grid=as.sparse, 
        transposed=transposed, .as.sparse=as.sparse, BPPARAM=BPPARAM)

    if (transposed) {
        rkout <- do.call(rbind, out)
    } else {
        rkout <- do.call(cbind, out)
    }

    if (withDimnames && !is.null(dimnames(x))) {
        dn <- list(rownames(x)[subset.row], colnames(x))
        if (transposed) { 
            dn <- rev(dn)
        }
        dimnames(rkout) <- dn
    } else if (!is.null(dimnames(rkout))) {
        dimnames(rkout) <- NULL
    }

    rkout
}

#' @importClassesFrom Matrix sparseMatrix dgCMatrix
#' @importFrom DelayedMatrixStats colRanks rowMins rowVars
#' @importFrom Matrix rowMeans
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.get_scaled_ranks <- function(block, transposed, .as.sparse) {
    if (is(block, "SparseArraySeed")) {
        block <- as(block, "sparseMatrix")
    }

    old <- getAutoBPPARAM()
    setAutoBPPARAM(NULL) # turning off any additional parallelization, just in case.
    on.exit(setAutoBPPARAM(old))

    out <- colRanks(DelayedArray(block), ties.method="average", preserveShape=FALSE)

    sig <- sqrt(rowVars(out) * (ncol(out)-1)) * 2
    if (any(is.na(sig) | sig==0)) {
        stop("rank variances of zero detected for a cell")
    }

    if (.as.sparse) {
        # Figure out what the zeroes got transformed into.
        is.zero <- which(block==0, arr.ind=TRUE)
        offset <- numeric(ncol(block))
        offset[is.zero[,2]] <- out[is.zero[,2:1]]

        out <- out - offset
        out <- as(out, "dgCMatrix")
    } else {
        out <- out - rowMeans(out)
    }

    out <- out/sig

    if (!transposed) {
        out <- t(out)
    }
    out
}
