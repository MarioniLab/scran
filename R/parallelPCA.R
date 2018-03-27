#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
.parallelPCA <- function(x, subset.row=NULL, scale=NULL, value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100,
                         niters=20, keep.perm=FALSE, approximate=FALSE, irlba.args=list(), BPPARAM=SerialParam())
# This performs Horn's parallel analysis to determine the number of PCs
# to retain, by randomizing each row and repeating the PCA to obtain
# an estimate of the mean variance explained per PC under a random model.
#
# written by Aaron Lun
# created 27 March 2018
{
    x0 <- x
    scale0 <- scale

    # Subsetting and scaling the matrix.
    if (!is.null(subset.row)) {
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        x <- x[subset.row,,drop=FALSE]
        scale <- scale[subset.row]
    }

    if (!is.null(scale)) {
        x <- x * scale
    }
    
    # Setting up the PCA function and its arguments.
    value <- match.arg(value)
    args <- list(y=t(x), max.rank=max.rank, value=value)

    if (approximate) {
        svdfun <- .irlba_svd
        args <- c(args, irlba.args)
    } else {
        svdfun <- .full_svd
    }

    # Running it once, and then multiple times after permutation.
    original <- do.call(svdfun, args)
    original.d2 <- original$d^2
    permuted <- bplapply(seq_len(niters), FUN=.parallel_PA, svdfun=svdfun, args=args, BPPARAM=BPPARAM)
    permuted.d2 <- Reduce("+", permuted)/niters

    # Figuring out where the original drops below the permuted.
    below <- permuted.d2 > original.d2
    if (!any(below)) {
        npcs <- max.rank
    } else {
        npcs <- min(which(below)) - 1L
    }
    npcs <- max(npcs, min.rank)

    # Actually computing the variance.
    var.exp <- original.d2 / (ncol(x) - 1)
    all.var <- sum(rowVars(DelayedArray(x)))
    
    # Figuring out what to return.
    out.val <- .convert_to_output(original, npcs, value, x0, scale0, subset.row)
    attr(out.val, "percentVar") <- var.exp/all.var
    
    if (keep.perm) { 
        permutations <- do.call(rbind, permuted.d2)/(ncol(x)-1L)
        permutations <- sweep(permutations, 2L, all.var, FUN="/")
        attr(out.val, "permuted.percentVar") <- permutations
    }

    return(out.val)
}

#' @importFrom stats prcomp
.full_svd <- function(y, max.rank, value) {
    y <- sweep(y, 2L, colMeans(y), check.margin = FALSE)

    if (value=="n") {
        nu <- nv <- 0
    } else if (value=="pcs") {
        nu <- max.rank
        nv <- 0
    } else {
        nu <- nv <- max.rank
    }

    out <- svd(y, nu=nu, nv=nv)
    out$d <- out$d[seq_len(max.rank)]
    return(out)
}

.irlba_svd <- function(y, max.rank, value, ...) {
    irlba::irlba(y, nv=max.rank, nu=max.rank, scale.=FALSE, center=TRUE, ...)
}

.parallel_PA <- function(svdfun, args, ...) {
    # Note that the ellipsis is ignored.
    args$y <- .Call(cxx_auto_shuffle, args$y)
    out <- do.call(svdfun, args)
    return(out$d^2)
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans
.convert_to_output <- function(svd.out, ncomp, value, x0, scale0, subset.row) {
    if (value=="n") {
        return(ncomp)
    } 
    
    ix <- seq_len(ncomp)
    U <- svd.out$u[,ix,drop=FALSE]
    D <- svd.out$d[ix]

    # Pulling out the PCs (column-multiplying the left eigenvectors).
    pcs <- sweep(U, 2L, D, FUN="*", check.margin = FALSE)
    if (value=="pcs") {
        rownames(pcs) <- colnames(x0)
        return(pcs)
    }

    # Creating a low-rank approximation.
    hits <- tcrossprod(pcs, V)

    if (!is.null(subset.row)) {
        output <- x0
        output[] <- 0
        output[subset.row,] <- hits

        # The idea is that after our SVD, we have X=UDV' where each column of X is a gene. 
        # Leftover genes are new columns in X, which are projected on the space of U by doing U'X.
        # This can be treated as new columns in DV', which can be multiplied by U to give denoised values.
        # I've done a lot of implicit transpositions here, hence the code does not tightly follow the logic above.
        leftovers <- !logical(nrow(x0))
        leftovers[subset.row] <- FALSE 
        left.x <- x0[leftovers,,drop=FALSE]
        left.scale <- scale0[leftovers]
        left.means <- rowMeans(DelayedArray(left.x))

        left.x <- left.x - left.mean
        left.x <- left.x * left.scale
        new.vals <- tcrossprod(left.x %*% U, U)
        new.vals <- new.vals / left.scale
        new.vals <- new.vals + left.means
        output[leftovers,] <- new.vals
    }

    return(output)
}
