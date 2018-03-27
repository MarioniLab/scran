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
        npcs <- length(below) 
    } else {
        npcs <- min(which(below)) - 1L
    }
    npcs <- .keep_rank_in_range(npcs, min.rank, length(original.d2))

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
.full_svd <- function(y, max.rank, value) 
# Convenience function for performing a SVD, with speed-ups
# to avoid computing the left and/or right eigenvectors if
# they are not necessary for the final 'value'.
{
    max.rank <- min(max.rank, dim(y))
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

.irlba_svd <- function(y, max.rank, value, extra.args=list()) 
# Convenience function for performing a SVD via the IRLBA,
# with protection against invalid inputs and increasing work
# to match the specified max.rank.
{
    arg.max <- pmatch(names(extra.args), "maxit")
    if (all(is.na(arg.max))) { 
        extra.args$maxit <- max(1000, max.rank*10)
    }
    max.rank <- min(max.rank, dim(y)-1L) # Note the -1 here, due to IRLBA's approximateness.
    all.args <- c(list(A=y, nv=max.rank, nu=max.rank, scale.=FALSE, center=TRUE), extra.args)
    do.call(irlba::irlba, all.args)
}

.parallel_PA <- function(svdfun, args, ...) 
# Function for use in bplapply, defined here to automatically take 
# advantage of the scran namespace when using snowParam.
{
    args$y <- .Call(cxx_auto_shuffle, args$y)
    out <- do.call(svdfun, args)
    return(out$d^2)
}

.keep_rank_in_range <- function(chosen, min.rank, nd)  
# A function to sensibly incorporate the min.rank information,
# while avoiding failures due to specification of a min.rank that is too large.
{
    max(chosen, min(min.rank, nd))
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans2
.convert_to_output <- function(svd.out, ncomp, value, original.mat, original.scale, subset.row) 
# Obtaining the desired output from the function; either the number of PCs,
# or the PCs themselves, or a low-rank approximation of the original matrix.
{
    if (value=="n") {
        return(ncomp)
    } 
    
    ix <- seq_len(ncomp)
    U <- svd.out$u[,ix,drop=FALSE]
    D <- svd.out$d[ix]

    # Pulling out the PCs (column-multiplying the left eigenvectors).
    pcs <- sweep(U, 2L, D, FUN="*", check.margin = FALSE)
    if (value=="pca") {
        colnames(pcs) <- sprintf("PC%i", ix)
        return(pcs)
    }

    # Creating a low-rank approximation by reforming UDV'. 
    # Note that we transpose to match dimensions, so it's really V(UD)'.
    V <- svd.out$v[,ix,drop=FALSE]
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
        if (!is.null(original.scale)) { 
            left.x <- left.x * original.scale[leftovers]
        }

        output[leftovers,] <- tcrossprod(left.x %*% U, U)
    }
    
    if (!is.null(original.scale)) {
        output <- output / original.scale
    }
    output <- output + all.means
    return(output)
}
