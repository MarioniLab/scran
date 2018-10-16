#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment assay
multiBatchPCA <- function(..., d=50, approximate=FALSE, irlba.args=list(), subset.row=NULL, assay.type="logcounts", get.spikes=FALSE, BPPARAM=SerialParam()) 
# Performs a multi-sample PCA (i.e., batches).
# Each batch is weighted inversely by the number of cells when computing the gene-gene covariance matrix.
# This avoids domination by samples with a large number of cells.
#
# written by Aaron Lun
# created 4 July 2018
{
    mat.list <- list(...)
    out <- .SCEs_to_matrices(mat.list, assay.type=assay.type, subset.row=subset.row, get.spikes=get.spikes)
    mat.list <- out$batches
    subset.row <- out$subset.row
    .multi_pca(mat.list, subset.row=subset.row, d=d, approximate=approximate, irlba.args=irlba.args, use.crossprod=TRUE, BPPARAM=BPPARAM) 
}

#' @importFrom DelayedArray DelayedArray t
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom BiocParallel SerialParam
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors metadata<-
.multi_pca <- function(mat.list, subset.row=NULL, d=50, approximate=FALSE, irlba.args=list(), use.crossprod=FALSE, BPPARAM=SerialParam()) 
# Internal function that uses DelayedArray to do the centering and scaling,
# to avoid actually realizing the matrices in memory.
{
    all.centers <- 0
    for (idx in seq_along(mat.list)) {
        current <- DelayedArray(mat.list[[idx]])
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }

        centers <- rowMeans2(current)
        all.centers <- all.centers + centers
        mat.list[[idx]] <- current
    }
    all.centers <- all.centers/length(mat.list) # grand average of centers (not using batch-specific centers, which makes compositional assumptions).

    centered <- scaled <- mat.list
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        current <- current - all.centers # centering each batch by the grand average.
        centered[[idx]] <- current
        current <- current/sqrt(ncol(current)) # downweighting samples with many cells.
        scaled[[idx]] <- t(current)
    }

    # Performing an SVD, if possible.
    if (d > min(ncol(scaled[[1]]), sum(vapply(scaled, FUN=nrow, FUN.VALUE=0L)))) {
        stop("'d' is too large for the number of cells and genes")
    }

    if (use.crossprod) {
        svd.out <- .fast_svd(scaled, nv=d, irlba.args=irlba.args, approximate=approximate, BPPARAM=BPPARAM)
    } else {
        combined <- as.matrix(do.call(rbind, scaled))
        if (!approximate) { 
            svd.out <- svd(combined, nu=0, nv=d)
        } else {
            svd.out <- do.call(irlba::irlba, c(list(A=combined, nu=0, nv=d), irlba.args))
        }
    }

    # Projecting the scaled matrices back into this space.
    final <- centered
    for (idx in seq_along(centered)) {
        final[[idx]] <- as.matrix(t(centered[[idx]]) %*% svd.out$v)
    }

    final <- as(final, "List")
    metadata(final) <- list(rotation=svd.out$v)
    return(final)
}

#' @importFrom DelayedArray t
#' @importFrom BiocParallel bplapply SerialParam
.fast_svd <- function(mat.list, nv, approximate=FALSE, irlba.args=list(), BPPARAM=SerialParam())
# Performs a quick irlba by performing the SVD on XtX or XXt,
# and then obtaining the V vector from one or the other.
{
    nrows <- sum(vapply(mat.list, FUN=nrow, FUN.VALUE=0L))
    ncols <- ncol(mat.list[[1]])
    
    # Creating the cross-product without actually rbinding the matrices.
    # This avoids creating a large temporary matrix.
    flipped <- nrows > ncols
    if (flipped) {
        collected <- bplapply(mat.list, FUN=.delayed_crossprod, BPPARAM=BPPARAM)
        final <- Reduce("+", collected)

    } else {
        final <- matrix(0, nrows, nrows)

        last1 <- 0L
        for (right in seq_along(mat.list)) {
            RHS <- mat.list[[right]]
            collected <- bplapply(mat.list[seq_len(right-1L)], FUN=.delayed_mult, Y=t(RHS), BPPARAM=BPPARAM)
            rdx <- last1 + seq_len(nrow(RHS))

            last2 <- 0L
            for (left in seq_along(collected)) {
                cross.prod <- collected[[left]]
                ldx <- last2 + seq_len(nrow(cross.prod))
                final[ldx,rdx] <- cross.prod
                final[rdx,ldx] <- t(cross.prod)
                last2 <- last2 + nrow(cross.prod)
            }
            
            last1 <- last1 + nrow(RHS)
        }

        tmat.list <- lapply(mat.list, t)
        diags <- bplapply(tmat.list, FUN=.delayed_crossprod, BPPARAM=BPPARAM)
        last1 <- 0L
        for (idx in seq_along(diags)) {
            indices <- last1 + seq_len(nrow(diags[[idx]]))
            final[indices,indices] <- diags[[idx]]
            last1 <- last1 + nrow(diags[[idx]])
        }
    }

    if (approximate) {
        svd.out <- do.call(irlba::irlba, c(list(A=final, nv=nv, nu=0), irlba.args))
    } else {
        svd.out <- svd(final, nv=nv, nu=0)
        svd.out$d <- svd.out$d[seq_len(nv)]
    }
    svd.out$d <- sqrt(svd.out$d)

    if (flipped) {
        # XtX means that the V in 'svd.out' is the original V,
        # which can be directly returned.
        return(svd.out)
    }

    # Otherwise, XXt means that the V in 'svd.out' is the original U.
    # We need to multiply the left-multiply the matrices by Ut, and then by D^-1.
    Ut <- t(svd.out$v)
    last <- 0L
    Vt <- matrix(0, nv, ncols)
    for (mdx in seq_along(mat.list)) {
        curmat <- mat.list[[mdx]]
        idx <- last + seq_len(nrow(curmat))
        Vt <- Vt + as.matrix(Ut[,idx,drop=FALSE] %*% curmat)
        last <- last + nrow(curmat)
    }
    Vt <- Vt / svd.out$d

    return(list(d=svd.out$d, v=t(Vt)))
}

.delayed_crossprod <- function(X, BPPARAM=SerialParam()) 
# DelayedMatrix crossprod, 1000 rows at a time.
{
    CHUNK <- 1000L
    last <- 0L
    output <- 0
    finish <- nrow(X)

    repeat {
        previous <- last + 1L
        last <- min(last + CHUNK, finish)
        block <- as.matrix(X[previous:last,,drop=FALSE])
        output <- output + crossprod(block)
        if (last==finish) { break }
    }

    return(output)
}

.delayed_mult <- function(X, Y, BPPARAM=SerialParam()) 
# DelayedMatrix multiplication, 1000 columns at a time.
{
    CHUNK <- 1000L
    last <- 0L
    output <- matrix(0, nrow(X), ncol(Y))
    finish <- ncol(Y)
    stopifnot(identical(ncol(X), nrow(Y)))

    repeat {
        previous <- min(last + 1L, finish)
        last <- min(last + CHUNK, finish)
        indices <- previous:last
        output[,indices] <- as.matrix(X %*% as.matrix(Y[,indices,drop=FALSE]))
        if (last==finish) { break }
    }

    return(output)
}
