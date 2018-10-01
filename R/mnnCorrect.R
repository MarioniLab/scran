#' @importFrom S4Vectors DataFrame Rle
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocGenerics t rbind
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans2
#' @export
mnnCorrect <- function(..., k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, 
                       svd.dim=0L, var.adj=TRUE, compute.angle=FALSE,
                       subset.row=NULL, order=NULL, pc.approx=FALSE, irlba.args=list(),
                       BPPARAM=SerialParam())
# Performs batch correction on multiple matrices of expression data,
# as specified in the ellipsis.
#    
# written by Laleh Haghverdi
# with modifications by Aaron Lun
# created 7 April 2017
{
    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }
    
    prep.out <- prepare.input.data(batches, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.out, subset.row=subset.row)
    in.batches <- prep.out$In
    out.batches <- prep.out$Out
    subset.row <- prep.out$Subset
    same.set <- prep.out$Same

    # Setting up the order.
    if (is.null(order)) {
        order <- seq_len(nbatches)
    } else {
        order <- as.integer(order)
        if (!identical(seq_len(nbatches), sort(order))) { 
            stop(sprintf("'order' should contain values in 1:%i", nbatches))
        }
    }
   
    # Setting up the variables.
    ref <- order[1]
    ref.batch.in <- t(in.batches[[ref]])
    if (!same.set) { 
        ref.batch.out <- t(out.batches[[ref]])
    }
    output <- vector("list", nbatches)
    output[[ref]] <- out.batches[[ref]]
    mnn.list <- vector("list", nbatches)
    mnn.list[[ref]] <- DataFrame(current.cell=integer(0), other.cell=integer(0), other.batch=integer(0))
    original.batch <- rep(ref, nrow(ref.batch.in)) 
    angle.list <- vector("list", nbatches)
    angle.list[[ref]] <- numeric(0)

    # Looping through the batches.
    for (b in 2:nbatches) { 
        target <- order[b]
        other.batch.in.untrans <- in.batches[[target]]
        other.batch.in <- t(other.batch.in.untrans)
        if (!same.set) { 
            other.batch.out.untrans <- out.batches[[target]]
            other.batch.out <- t(other.batch.out.untrans)
        }
        
        # Finding pairs of mutual nearest neighbours.
        sets <- find.mutual.nn(ref.batch.in, other.batch.in, k1=k, k2=k, BPPARAM=BPPARAM)
        s1 <- sets$first
        s2 <- sets$second      

        # Computing the correction vector.
        correction.in <- compute.correction.vectors(ref.batch.in, other.batch.in, s1, s2, other.batch.in.untrans, sigma)
        if (!same.set) {
            correction.out <- compute.correction.vectors(ref.batch.out, other.batch.out, s1, s2, other.batch.in.untrans, sigma)
            # NOTE: use of 'other.batch.in.untrans' here is intentional, 
            # as the distances should be the same as the MNN distances.
        }

        # Calculating the smallest angle between each correction vector and the first 2 basis vectors of the reference batch.
        if (compute.angle) {
            ref.centred <- t(ref.batch.in)
            ref.centred <- ref.centred - rowMeans2(DelayedArray(ref.centred))
            ref.basis <- .svd_internal(ref.centred, nu=2, nv=0, pc.approx=pc.approx, irlba.args=irlba.args)$u

            angle.out <- numeric(nrow(correction.in))
            for (i in seq_along(angle.out)) {
                angle.out[i] <- find.shared.subspace(ref.basis, t(correction.in[i,,drop=FALSE]))$angle
            }
            angle.list[[target]] <- angle.out
        }

        # Removing any component of the correction vector that's parallel to the biological basis vectors in either batch.
        # Only using cells involved in MNN pairs, to avoid undercorrection in directions that weren't problematic anyway
        # (given that using SVDs was intended to mitigate the effect of identifying the wrong MNN pairs).
        if (svd.dim>0) {
            u1 <- unique(s1)
            u2 <- unique(s2)

            # Computing the biological subspace in both batches, and subtract it from the batch correction vector.
            in.span1 <- get.bio.span(t(ref.batch.in[u1,,drop=FALSE]), ndim=svd.dim,
                                     pc.approx=pc.approx, irlba.args=irlba.args)
            in.span2 <- get.bio.span(other.batch.in.untrans[,u2,drop=FALSE], ndim=svd.dim, 
                                     pc.approx=pc.approx, irlba.args=irlba.args)
            correction.in <- subtract.bio(correction.in, in.span1, in.span2)

            # Repeating for the output values.
            if (!same.set) { 
                out.span1 <- get.bio.span(t(ref.batch.out[u1,,drop=FALSE]), subset.row=subset.row,
                                          ndim=svd.dim, pc.approx=pc.approx, irlba.args=irlba.args)
                out.span2 <- get.bio.span(other.batch.out.untrans[,u2,drop=FALSE], subset.row=subset.row,
                                          ndim=svd.dim, pc.approx=pc.approx, irlba.args=irlba.args)
                correction.out <- subtract.bio(correction.out, out.span1, out.span2, subset.row=subset.row)
            }
        } 
       
        # Adjusting the shift variance; done after any SVD so that variance along the correction vector is purely technical.
        if (var.adj) { 
            correction.in <- adjust.shift.variance(ref.batch.in, other.batch.in, correction.in, sigma=sigma)
            if (!same.set) {
                correction.out <- adjust.shift.variance(ref.batch.out, other.batch.out, correction.out, sigma=sigma, subset.row=subset.row) 
            }
        }

        # Applying the correction and expanding the reference batch. 
        other.batch.in <- other.batch.in + correction.in
        ref.batch.in <- rbind(ref.batch.in, other.batch.in)
        if (same.set) {
            output[[target]] <- t(other.batch.in)
        } else {
            other.batch.out <- other.batch.out + correction.out
            ref.batch.out <- rbind(ref.batch.out, other.batch.out)
            output[[target]] <- t(other.batch.out)
        }

        # Storing the identities of the MNN pairs (using RLEs for compression of runs).
        mnn.list[[target]] <- DataFrame(current.cell=s2, other.cell=Rle(s1), other.batch=Rle(original.batch[s1]))
        original.batch <- c(original.batch, rep(target, nrow(other.batch.in)))
    }

    # Formatting output to be consistent with input.
    names(output) <- names(batches)
    names(mnn.list) <- names(batches)
    final <- list(corrected=output, pairs=mnn.list)
    if (compute.angle){ 
        names(angle.list) <- names(batches)
        final$angles <- angle.list
    }
    return(final)
}

prepare.input.data <- function(batches, cos.norm.in, cos.norm.out, subset.row) {
    nbatches <- length(batches)

    # Checking for identical number of rows (and rownames).
    first <- batches[[1]]
    ref.nrow <- nrow(first)
    ref.rownames <- rownames(first)
    for (b in 2:nbatches) {
        current <- batches[[b]]
        if (!identical(nrow(current), ref.nrow)) {
            stop("number of rows is not the same across batches")
        } else if (!identical(rownames(current), ref.rownames)) {
            stop("row names are not the same across batches")
        }
    }

    # Subsetting to the desired subset of genes.
    in.batches <- out.batches <- batches
    same.set <- TRUE
    if (!is.null(subset.row)) { 
        subset.row <- .subset_to_index(subset.row, batches[[1]], byrow=TRUE)
        if (identical(subset.row, seq_len(ref.nrow))) { 
            subset.row <- NULL
        } else {
            same.set <- FALSE
            in.batches <- lapply(in.batches, "[", i=subset.row, , drop=FALSE) # Need the extra comma!
        }
    }

    # Applying cosine normalization for MNN identification. 
    # We use the L2 norm for the subsetted input to adjust the output, 
    # to ensure that results are consistent regardless of the manner of subsetting.
    if (cos.norm.in) { 
        norm.scaling <- vector("list", nbatches)
        for (b in seq_len(nbatches)) { 
            current.in <- in.batches[[b]]
            cos.out <- cosineNorm(current.in, mode="all")
            in.batches[[b]] <- cos.out$matrix
            norm.scaling[[b]] <- cos.out$l2norm
        }
    }
    if (cos.norm.out) { 
        if (!cos.norm.in) { 
            norm.scaling <- lapply(in.batches, cosineNorm, mode="l2norm")
        }
        for (b in seq_len(nbatches)) { 
            out.batches[[b]] <- sweep(out.batches[[b]], 2, norm.scaling[[b]], "/")
        }
    }
    if (cos.norm.out!=cos.norm.in) { 
        same.set <- FALSE
    }

    return(list(In=in.batches, Out=out.batches, Subset=subset.row, Same=same.set))
}

#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
find.mutual.nn <- function(data1, data2, k1, k2, BNPARAM=NULL, BPPARAM=SerialParam()) 
# Finds mutal neighbors between data1 and data2.
{
    data1 <- as.matrix(data1)
    data2 <- as.matrix(data2)
    W21 <- queryKNN(data2, query=data1, k=k1, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=FALSE)
    W12 <- queryKNN(data1, query=data2, k=k2, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=FALSE)
    out <- .Call(cxx_find_mutual_nns, W21$index, W12$index)
    names(out) <- c("first", "second")
    return(out)
}

compute.correction.vectors <- function(data1, data2, mnn1, mnn2, original2, sigma) 
# Computes the batch correction vector for each cell in 'data2'.
# 'original2' should also be supplied to compute distances 
# (this may not be the same as 't(data2)' due to normalization, subsetting).
{      
    vect <- data1[mnn1,,drop=FALSE] - data2[mnn2,,drop=FALSE]    
    cell.vect <- .Call(cxx_smooth_gaussian_kernel, vect, mnn2-1L, original2, sigma)
    return(t(cell.vect)) 
}

adjust.shift.variance <- function(data1, data2, correction, sigma, subset.row=NULL) 
# Performs variance adjustment to avoid kissing effects.    
{
    cell.vect <- correction 
    if (!is.null(subset.row)) { 
        # Only using subsetted genes to compute locations, consistent with our policy in SVD. 
        cell.vect <- cell.vect[,subset.row,drop=FALSE]
        data1 <- data1[,subset.row,drop=FALSE]
        data2 <- data2[,subset.row,drop=FALSE]
    }
    scaling <- .Call(cxx_adjust_shift_variance, t(data1), t(data2), cell.vect, sigma)
    scaling <- pmax(scaling, 1)
    return(scaling * correction)
}

.svd_internal <- function(X, nu, nv, pc.approx=FALSE, irlba.args=list()) {
    if (!pc.approx) { 
        S <- svd(X, nu=nu, nv=nv)
    } else {
        # Note: can't use centers=, as that subtracts from columns not rows.
        S <- do.call(irlba::irlba, c(list(A=X, nu=nu, nv=nv), irlba.args))
    }
    return(S)
}

get.bio.span <- function(exprs, ndim, subset.row=NULL, pc.approx=FALSE, irlba.args=list())
# Computes the basis matrix of the biological subspace of 'exprs'.
# The first 'ndim' dimensions are assumed to capture the biological subspace.
{
    centred <- exprs - rowMeans(exprs)
    centred <- as.matrix(centred)

    if (!is.null(subset.row)) {
        leftovers <- centred[-subset.row,,drop=FALSE]
        centred <- centred[subset.row,,drop=FALSE]
        nv <- ndim
    } else {
        nv <- 0
    }

    # Returning the basis vectors directly if there was no subsetting. 
    ndim <- min(ndim, dim(centred))
    nv <- min(nv, ndim)
    S <- .svd_internal(centred, nu=ndim, nv=nv, pc.approx=pc.approx, irlba.args=irlba.args)
    if (is.null(subset.row)) { 
        return(S$u)
    } 

    # Computing the basis values for the unused genes.
    output <- matrix(0, nrow(exprs), ndim)
    output[subset.row,] <- S$u
    leftovers <- leftovers %*% S$v 
    leftovers <- sweep(leftovers, 2, S$d[seq_len(ndim)], "/")
    output[-subset.row,] <- leftovers
    return(output)
}

subtract.bio <- function(correction, span1, span2, subset.row=NULL) 
# Computes the component parallel biological basis vectors in the correction vectors,
# and subtracts them. Note that the order of span1 and span2 does not matter.
{ 
    for (span in list(span1, span2)) { 
        if (is.null(subset.row)) { 
            bio.mag <- correction %*% span
        } else { 
            # Only using the subset to compute the magnitude that is parallel.
            bio.mag <- correction[,subset.row,drop=FALSE] %*% span[subset.row,,drop=FALSE]
        }
        bio.comp <- bio.mag %*% t(span)
        correction <- correction - bio.comp
    }
    return(correction)
}    

find.shared.subspace <- function(A, B, sin.threshold=0.85, cos.threshold=1/sqrt(2), 
                                 assume.orthonormal=FALSE, get.angle=TRUE) 
# Computes the maximum angle between subspaces, to determine if spaces are orthogonal.
# Also identifies if there are any shared subspaces. 
{
    if (!assume.orthonormal) { 
        A <- pracma::orth(A)
        B <- pracma::orth(B)
    }
     
    # Singular values close to 1 indicate shared subspace A \invertsymbol{U} B
    # Otherwise A and B are completely orthogonal, i.e., all angles=90.
    cp.AB <- crossprod(A, B)
    S <- svd(cp.AB)
    shared <- sum(S$d > sin.threshold)
    if (!get.angle) {
        return(list(nshared=shared))
    }

    # Computing the angle from singular values; using cosine for large angles,
    # sine for small angles (due to differences in relative accuracy).
    costheta <- min(S$d) 
    if (costheta < cos.threshold){ 
        theta <- acos(min(1, costheta))
    } else {
        if (ncol(A) < ncol(B)){ 
            sintheta <- svd(t(A) - tcrossprod(cp.AB, B), nu=0, nv=0)$d[1]
        } else {
            sintheta <- svd(t(B) - tcrossprod(t(cp.AB), A), nu=0, nv=0)$d[1]
        }
        theta <- asin(min(1, sintheta)) 
    }
    
    list(angle=180*theta/pi, nshared=shared)
}
