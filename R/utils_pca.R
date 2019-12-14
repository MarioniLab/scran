#' @importFrom BiocSingular runSVD bsparam
#' @importFrom BiocParallel SerialParam
.centered_SVD <- function(y, max.rank, BSPARAM=bsparam(), BPPARAM=SerialParam(), keep.left=TRUE, keep.right=TRUE)
# Performs the PCA given a log-expression matrix.
# Switches between svd() and irlba() on request.
# Output format is guaranteed to be the same.
{
    runSVD(y, center=TRUE, BSPARAM=BSPARAM, k=max.rank, 
        nu=if (keep.left) max.rank else 0L,
        nv=if (keep.right) max.rank else 0L,
        BPPARAM=BPPARAM)
}

.keep_rank_in_range <- function(chosen, min.rank, nd)
# A function to sensibly incorporate the min.rank information,
# while avoiding failures due to specification of a min.rank that is too large.
{
    max(chosen, min(min.rank, nd))
}

.svd_to_pca <- function(svd.out, ncomp, named=TRUE) 
# Converts centred results to PCs.
{
    if (is.null(svd.out$u)) {
        stop("missing 'U' in SVD results")
    } else if (ncomp > ncol(svd.out$u)) {
        warning("requested number of components greater than available rank")
        ncomp <- ncol(svd.out$u)
    }

    ix <- seq_len(ncomp)
    U <- svd.out$u[,ix,drop=FALSE]
    D <- svd.out$d[ix]

    # Pulling out the PCs (column-multiplying the left eigenvectors).
    pcs <- sweep(U, 2L, D, FUN="*", check.margin = FALSE)
    if (named) {
        colnames(pcs) <- sprintf("PC%i", ix)
    }
    return(pcs)
}
