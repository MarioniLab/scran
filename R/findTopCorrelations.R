#' @importFrom Matrix t
#' @importFrom BiocSingular IrlbaParam
#' @importFrom BiocNeighbors findKNN queryKNN KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom stats pt p.adjust
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray is_sparse
#' @importFrom beachmat rowBlockApply
.find_top_correlations <- function(x, k=10, y=NULL, d=50, use.names=TRUE, BSPARAM=IrlbaParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) {
    if (is.null(y)) {
        names1 <- names2 <- rownames(x)
        rank.x <- .create_rank_matrix(t(x), deferred=is_sparse(x), BPPARAM=BPPARAM, transposed=TRUE)

        search.x <- .compress_rank_matrix(rank.x, d=d, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        nn.out <- findKNN(search.x, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

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
            k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

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
    if (!is.null(d)) {
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
