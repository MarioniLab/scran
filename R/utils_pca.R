#' @importFrom BiocSingular runSVD IrlbaParam ExactParam
#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix colMeans
.centered_SVD <- function(y, max.rank, BSPARAM=ExactParam(), BPPARAM=SerialParam(), approximate=FALSE, extra.args=list(), keep.left=TRUE, keep.right=TRUE)
# Performs the PCA given a log-expression matrix.
# Switches between svd() and irlba() on request.
# Output format is guaranteed to be the same.
{
    if (!is.null(approximate)) {
        .Deprecated(msg="'approximate=TRUE' or 'pc.approx=TRUE' are deprecated.\nUse BSPARAM=BiocSingular::IrlbaParam() instead.")
        if (approximate) {
            BSPARAM <- do.call(IrlbaParam, extra.args)
        } else {
            BSPARAM <- ExactParam()
        }
    }

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

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans2
.svd_to_lowrank <- function(svd.out, ncomp, original.mat, subset.row)
# Obtaining the desired output from the function; either the number of PCs,
# or the PCs themselves, or a low-rank approximation of the original matrix.
{
    if (is.null(svd.out$u) || is.null(svd.out$v)) {
        stop("missing 'U' or 'V' in SVD results")
    } else if (ncomp > ncol(svd.out$u)) {
        warning("requested number of components greater than available rank")
        ncomp <- ncol(svd.out$u)
    }

    ix <- seq_len(ncomp)
    U <- svd.out$u[,ix,drop=FALSE]
    V <- svd.out$v[,ix,drop=FALSE]
    D <- svd.out$d[ix]

    # Creating a low-rank approximation by reforming UDV'.
    # Note that we transpose to match dimensions, so it's really V(UD)'.
    pcs <- sweep(U, 2L, D, FUN="*", check.margin = FALSE)
    hits <- tcrossprod(V, pcs)
    all.means <- rowMeans2(DelayedArray(original.mat))

    if (is.null(subset.row)) {
        output <- hits
        dimnames(output) <- dimnames(original.mat)
    } else {
        output <- original.mat
        output[subset.row,] <- hits

        # The idea is that after our SVD, we have X=UDV' where each column of X is a gene.
        # Leftover genes are new columns in X, which are projected on the space of U by doing U'X.
        # This can be treated as new columns in DV', which can be multiplied by U to give denoised values.
        # I've done a lot of implicit transpositions here, hence the code does not tightly follow the logic above.
        leftovers <- !logical(nrow(original.mat))
        leftovers[subset.row] <- FALSE

        left.x <- original.mat[leftovers,,drop=FALSE] - all.means[leftovers]

        output[leftovers,] <- tcrossprod(left.x %*% U, U)
    }

    return(output + all.means)
}
