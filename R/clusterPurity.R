#' Evaluate cluster purity
#'
#' 
NULL

#' @importFrom BiocNeighbors KmknnParam buildIndex findKNN findNeighbors
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scater .bpNotSharedOrUp
#' @importFrom Matrix t
.cluster_purity <- function(x, clusters, k=50, transposed=FALSE, subset.row=NULL,
    pseudo.count=1L, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    if (!transposed) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        x <- t(x)
    }

    x <- as.matrix(x)
    idx <- buildIndex(x, BNPARAM=BNPARAM)

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    dist <- findKNN(BNINDEX=idx, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, last=1, get.index=FALSE)
    nout <- findNeighbors(BNINDEX=idx, threshold=dist, BNPARAM=BNPARAM, get.distance=FALSE)$index

    # Counting the number of entities in the same cluster as me.
    sameness <- clusters[unlist(nout)]==rep(clusters, lengths(nout))
    samelist <- relist(sameness, nout)
    nsame <- vapply(samelist, sum, 0L)

    (nsame + pseudo.count)/(lengths(samelist) + 2*pseudo.count)
}
