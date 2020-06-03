#' Cell cycle phase training
#'
#' Use gene expression data to train a classifier for cell cycle phase.
#' 
#' @param x A numeric matrix of gene expression values where rows are genes and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param phases A list of subsetting vectors specifying which cells are in each phase of the cell cycle.
#' This should typically be of length 3, with elements named as \code{"G1"}, \code{"S"} and \code{"G2M"}.
#' @param gene.names A character vector of gene names.
#' @param fraction A numeric scalar specifying the minimum fraction to define a marker gene pair.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param ... For the generic, additional arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, additional arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.
#' 
#' @details
#' This function implements the training step of the pair-based prediction method described by Scialdone et al. (2015).
#' Pairs of genes (A, B) are identified from a training data set where in each pair,
#'     the fraction of cells in phase G1 with expression of A > B (based on expression values in \code{training.data}) 
#'     and the fraction with B > A in each other phase exceeds \code{fraction}.
#' These pairs are defined as the marker pairs for G1.
#' This is repeated for each phase to obtain a separate marker pair set.
#' 
#' Pre-defined sets of marker pairs are provided for mouse and human (see Examples).
#' The mouse set was generated as described by Scialdone et al. (2015), while the human training set was generated with data from Leng et al. (2015).
#' Classification from test data can be performed using the \code{\link{cyclone}} function.
#' For each cell, this involves comparing expression values between genes in each marker pair. 
#' The cell is then assigned to the phase that is consistent with the direction of the difference in expression in the majority of pairs.
#' 
#' @return
#' A named list of data.frames, where each data frame corresponds to a cell cycle phase and contains the names of the genes in each marker pair.
#' 
#' @author
#' Antonio Scialdone,
#' with modifications by Aaron Lun
#' 
#' @seealso
#' \code{\link{cyclone}}, to perform the classification on a test dataset.
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ncells=50, ngenes=200)
#' 
#' is.G1 <- 1:20
#' is.S <- 21:30
#' is.G2M <- 31:50
#' out <- sandbag(sce, list(G1=is.G1, S=is.S, G2M=is.G2M))
#' str(out)
#' 
#' # Getting pre-trained marker sets
#' mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
#' hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
#' 
#' @references
#' Scialdone A, Natarajana KN, Saraiva LR et al. (2015). 
#' Computational assignment of cell-cycle stage from single-cell transcriptome data.
#' \emph{Methods} 85:54--61
#' 
#' Leng N, Chu LF, Barry C et al. (2015).
#' Oscope identifies oscillatory genes in unsynchronized single-cell RNA-seq experiments.
#' \emph{Nat. Methods} 12:947--50
#' 
#' @name sandbag
NULL

find.markers <- function(current.data, other.data, gene.names, fraction=0.5)
# This identifies pairs of genes whose relative expression is > 0 in 
# at least a 'fraction' of cells in one phase is < 0 in at least 
# 'fraction' of the cells in each of the other phases.
{
    Ngenes <- ncol(current.data)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' vector must be equal to 'x' nrows")
    }
    other.Ngenes <- vapply(other.data, ncol, FUN.VALUE=0L)
    if (any(other.Ngenes!=Ngenes)) { 
        stop("number of genes in each class must be the same")
    }
    Ncells <- nrow(current.data)
    other.Ncells <- vapply(other.data, nrow, FUN.VALUE=0L)
    if (Ncells==0L || any(other.Ncells==0L)) {
        stop("each class must have at least one cell")
    }

    # Calculating thresholds.
    Nthr.cur <- ceiling(Ncells * fraction)
    Nthr.other <- ceiling(other.Ncells * fraction)

    if (Ngenes) { 
        collected <- vector("list", Ngenes*2)
        collected[[1]] <- matrix(0L, 0, 2)

        for (i in seq_len(Ngenes-1L)) { 
            others <- (i+1):Ngenes
            cur.diff <- current.data[,i] - current.data[,others,drop=FALSE]
            other.diff <- lapply(other.data, function(odata) { odata[,i] - odata[,others,drop=FALSE] })

            # Looking for marker pairs that are up in the current group and down in the other groups.
            cur.pos.above.threshold <- colSums(cur.diff > 0) >= Nthr.cur
            other.neg.above.threshold <- mapply(function(odata, thr) { colSums(odata < 0) >= thr}, 
                                                other.diff, Nthr.other, SIMPLIFY=FALSE, USE.NAMES=FALSE)
            chosen <- others[cur.pos.above.threshold & Reduce(`&`, other.neg.above.threshold)]
            if (length(chosen)) { 
                collected[[i*2]] <- cbind(i, chosen)
            }

            # Looking for marker pairs that are down in the current group and up in the other groups.
            cur.neg.above.threshold <- colSums(cur.diff < 0) >= Nthr.cur
            other.pos.above.threshold <- mapply(function(odata, thr) { colSums(odata > 0) >= thr}, 
                                                other.diff, Nthr.other, SIMPLIFY=FALSE, USE.NAMES=FALSE)
            chosen.flip <- others[cur.neg.above.threshold & Reduce(`&`, other.pos.above.threshold)]
            if (length(chosen.flip)) { 
                collected[[i*2+1]] <- cbind(chosen.flip, i)
            }
        }

        collected <- do.call(rbind, collected)
        g1 <- gene.names[collected[,1]]
        g2 <- gene.names[collected[,2]]
    } else {
        g1 <- g2 <- character(0)
    }

    data.frame(first=g1, second=g2, stringsAsFactors=FALSE)
}

#' @export
#' @rdname sandbag
setGeneric("sandbag", function(x, ...) standardGeneric("sandbag"))

#' @export
#' @rdname sandbag
#' @importFrom scuttle .subset2index
setMethod("sandbag", "ANY", function(x, phases, gene.names=rownames(x), fraction=0.5, subset.row=NULL) {
    subset.row <- .subset2index(subset.row, x, byrow=TRUE)
    gene.names <- gene.names[subset.row]

    class.names <- names(phases)
    if (is.null(class.names) || any(is.na(class.names))) {
        stop("'phases' must have non-missing, non-NULL names") 
    }

    gene.data <- lapply(phases, function(cl) t(x[subset.row,cl,drop=FALSE]))
    nclasses <- length(gene.data)
    marker.pairs <- vector("list", nclasses)
    for (i in seq_len(nclasses)) {
        marker.pairs[[i]] <- find.markers(gene.data[[i]], gene.data[-i], fraction=fraction, gene.names=gene.names)
    }

    names(marker.pairs) <- class.names
    marker.pairs
})

#' @export
#' @rdname sandbag
#' @importFrom SummarizedExperiment assay
setMethod("sandbag", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    sandbag(assay(x, i=assay.type), ...)
})
