#' Find top correlations between genes
#'
#' For each gene, find the subset of other genes that have strongest positive/negative Spearman's rank correlations in a normalized expression matrix.
#'
#' @param x A numeric matrix containing a (possibly log-transformed) normalized expression matrix for genes (rows) and cells (columns).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param number Integer scalar specifying the number of top correlated genes to report for each gene in \code{x}.
#' @param y Another normalized expression matrix with the same number of cells as \code{x}, possibly with different features.
#' For the SummarizedExperiment method, this may also be a SummarizedExperiment object containing such a matrix.
#' @param d Integer scalar specifying the number of dimensions to use for the approximate search via PCA.
#' If \code{NA}, no approximation of the rank values is performed prior to the search.
#' @param use.names Logical scalar specifying whether row names of \code{x} and/or \code{y} should be reported in the output, if available.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for the PCA.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the neighbor search.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization scheme to use.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.type String or integer scalar specifying the assay containing the matrix of interest in \code{x} (and \code{y}, if a SummarizedExperiment).
#'
#' @return A \linkS4class{DataFrame} with up to \code{nrow(x) * number} rows, containing the top \code{number} correlated genes for each gene in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{gene1}, the name (character) or row index (integer) of each gene in \code{x}.
#' \item \code{gene2}, the name (character) or row index (integer) of one of the top correlated genes to \code{gene1}.
#' This is another gene in \code{x} if \code{y=NULL}, otherwise it is a gene in \code{y}.
#' \item \code{rho}, the Spearman rank correlation for the current pair of \code{gene1} and \code{gene2}.
#' \item \code{p.value}, the approximate p-value associated with \code{rho} under the null hypothesis that the correlation is zero.
#' \item \code{FDR}, the adjusted p-value.
#' }
#' 
#' @details
#' In most cases, we only care about the top-correlated genes, allowing us to skip a lot of unnecessary computation.
#' This is achieved by transforming the problem of finding the largest Spearman correlation into a nearest-neighbor search in rank space.
#' For the sake of speed, we approximate the search by performing PCA to compress the rank values for all genes.
#'
#' We compute the p-value for each gene using the approximate method implemented in \code{\link{cor.test}}.
#' The FDR correction is performed by considering all possible pairs of genes, as these are implicitly tested in the neighbor search.
#' Note that this is somewhat conservative as it does not consider strong correlations outside the reported genes.
#'
#' @author Aaron Lun
#'
#' @name findTopCorrelations
NULL

######################
##### Internals ######
######################

#' @importFrom Matrix t
#' @importFrom BiocSingular IrlbaParam
#' @importFrom BiocNeighbors findKNN queryKNN KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pt p.adjust
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray is_sparse
#' @importFrom beachmat rowBlockApply
.find_top_correlations <- function(x, number=10, y=NULL, d=50, use.names=TRUE, BSPARAM=IrlbaParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) {
    if (is.null(y)) {
        names1 <- names2 <- rownames(x)
        rank.x <- .create_rank_matrix(t(x), deferred=is_sparse(x), BPPARAM=BPPARAM, transposed=TRUE)

        search.x <- .compress_rank_matrix(rank.x, d=d, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        nn.out <- findKNN(search.x, k=number, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

        rho <- rowBlockApply(rank.x, FUN=.compute_exact_neighbor_distances, 
            other=rank.x, indices=nn.out$index, BPPARAM=BPPARAM)
        rho <- do.call(rbind, rho)

        ntests <- nrow(x) * (nrow(x) - 1L) / 2L

    } else {
        names1 <- rownames(x)
        names2 <- rownames(y)
        combined <- rbind(t(x), t(y))
        rank.com <- .create_rank_matrix(combined, deferred=is_sparse(combined), BPPARAM=BPPARAM, transposed=TRUE)

        search.com <- .compress_rank_matrix(rank.x, d=d, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        in.first <- seq_len(nrow(x))
        in.second <- nrow(x) + seq_len(nrow(y))
        nn.out <- queryKNN(X=search.com[in.second,,drop=FALSE], query=search.com[in.first,,drop=FALSE],
            k=number, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

        rho <- rowBlockApply(rank.com[in.first,,drop=FALSE], FUN=.compute_exact_neighbor_distances, 
            other=rank.com[in.second,,drop=FALSE], indices=nn.out$index, BPPARAM=BPPARAM)
        rho <- do.call(rbind, rho)

        ntests <- nrow(x) * nrow(y)
    }

    self <- rep(seq_len(nrow(x)), ncol(nn.out$index))
    df <- DataFrame(gene1=self, gene2=as.vector(nn.out$index))

    # Assembling the DF.
    if (use.names) {
        if (!is.null(names2)) {
            df$gene2 <- names2[df$gene2]
        } 
        if (!is.null(names1)) {
            df$gene1 <- names1[df$gene1]
        }
    }

    # Mildly adapted from cor.test.
    ncells <- ncol(x)
    r <- 1 - df$rho/(ncells * (ncells^2 - 1)/6)
    df$p.value <- pt(r/sqrt((1 - r^2)/(ncells - 2)), df = ncells - 2, lower.tail = (df$rho <= (ncells^3 - ncells)/6))
    df$FDR <- p.adjust(df$p.value, method="BH", n = ntests)

    o <- order(df$self, df$p.value)
    df[o,,drop=FALSE]
}

.compress_rank_matrix <- function(rank.x, d, BSPARAM, BPPARAM) {
    if (!is.na(d)) {
        svd.out <- .centered_SVD(rank.x, max.rank=d, keep.right=FALSE, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        .svd_to_pca(svd.out, d, named=FALSE)
    } else {
        as.matrix(rank.x)
    }
}

#' @importFrom DelayedArray currentViewport makeNindexFromArrayViewport
.compute_exact_neighbor_distances <- function(block, other, indices) {
    vp <- currentViewport()
    subset <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)
    if (!is.null(subset[[1]])) {
        indices <- indices[subset[[1]],,drop=FALSE]
    }

    output <- matrix(0, nrow(block), ncol(indices))
    for (j in seq_len(ncol(indices))) {
        current <- block - as.matrix(rank.x[indices[,j],,drop=FALSE])
        output[,j] <- sqrt(rowSums(current^2))
    }

    output
}

#######################
##### S4 methods ######
#######################

#' @export
#' @rdname findTopCorrelations
setGeneric("findTopCorrelations", function(x, number, ...) standardGeneric("findTopCorrelations"))

#' @export
#' @rdname findTopCorrelations
setMethod("findTopCorrelations", "ANY", .find_top_correlations)

#' @export
#' @rdname findTopCorrelations
#' @importFrom SummarizedExperiment assay
setMethod("findTopCorrelations", "SummarizedExperiment", function(x, number, y=NULL, ..., assay.type="logcounts") {
    if (is(y, "SummarizedExperiment")) {
        y <- assay(y, assay.type)
    }
    x <- assay(x, assay.type)
    .find_top_correlations(x, number=number, y=y, ...)
})
