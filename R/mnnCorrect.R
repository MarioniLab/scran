mnnCorrect <- function(..., k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, svd.dim=0L, var.adj=TRUE,
                       subset.row=NULL, order=NULL, pc.approx=FALSE, BPPARAM=SerialParam())
# Performs batch correction on multiple matrices of expression data,
# as specified in the ellipsis.
#    
# written by Laleh Haghverdi
# with modifications by Aaron Lun
# created 7 April 2017
# last modified 6 November 2017
{
    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }
   
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
            cos.out <- cosine.norm(current.in, mode="all")
            in.batches[[b]] <- cos.out$matrix
            norm.scaling[[b]] <- cos.out$l2norm
        }
    }
    if (cos.norm.out) { 
        if (!cos.norm.in) { 
            norm.scaling <- lapply(in.batches, cosine.norm, mode="l2norm")
        }
        for (b in seq_len(nbatches)) { 
            out.batches[[b]] <- sweep(out.batches[[b]], 2, norm.scaling[[b]], "/")
        }
    }
    if (cos.norm.out!=cos.norm.in) { 
        same.set <- FALSE
    }

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

        if (svd.dim>0) {
            u1 <- unique(s1)
            u2 <- unique(s2)

            # Computing the biological subspace in both batches.
            in.ndim <- min(c(svd.dim, dim(ref.batch.in), dim(other.batch.in)))
            in.span1 <- get.bio.span(t(ref.batch.in[u1,,drop=FALSE]), ndim=min(in.ndim, length(s1)), pc.approx=pc.approx)
            in.span2 <- get.bio.span(other.batch.in.untrans[,u2,drop=FALSE], ndim=min(in.ndim, length(s2)), pc.approx=pc.approx)

            # Reduce the component in each span from the batch correction vector.
            correction.in <- subtract.bio(correction.in, in.span1, in.span2)

            # Repeating for the output values.
            if (!same.set) { 
                out.ndim <- min(c(svd.dim, dim(ref.batch.out), dim(other.batch.out)))
                if (!is.null(subset.row)) { 
                    out.ndim <- min(out.ndim, length(subset.row))
                }                    
                
                out.span1 <- get.bio.span(t(ref.batch.out[u1,,drop=FALSE]), subset.row=subset.row,
                                          ndim=min(out.ndim, length(s1)), pc.approx=pc.approx)
                out.span2 <- get.bio.span(other.batch.out.untrans[,u2,drop=FALSE], subset.row=subset.row,
                                          ndim=min(out.ndim, length(s2)), pc.approx=pc.approx)
                correction.out <- subtract.bio(correction.out, out.span1, out.span2, subset.row=subset.row)
            }
        } 
       
        # Adjusting the shift variance, if requested.
        # This is done after any SVD so that variance along the correction vector is purely technical.
        if (var.adj) { 
            correction.in <- adjust.shift.variance(ref.batch.in, other.batch.in, correction.in)
            if (!same.set) {
                correction.out <- adjust.shift.variance(ref.batch.out, other.batch.out, correction.out, subset.row=subset.row) 
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
    list(corrected=output, pairs=mnn.list)
}

find.mutual.nn <- function(data1, data2, k1, k2, BPPARAM) 
# Finds mutal neighbors between data1 and data2.
{
    W21 <- bpl.get.knnx(data2, query=data1, k=k1, BPPARAM=BPPARAM)
    W12 <- bpl.get.knnx(data1, query=data2, k=k2, BPPARAM=BPPARAM)
    out <- .Call(cxx_find_mutual_nns, W21$nn.index, W12$nn.index)
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

adjust.shift.variance <- function(data1, data2, correction, subset.row=NULL) 
# Performs variance adjustment to avoid kissing effects.    
{
    # Only using subsetted genes to compute locations, consistent with our policy in SVD. 
    cell.vect <- correction 
    if (!is.null(subset.row)) { 
        cell.vect <- cell.vect[,subset.row,drop=FALSE]
        data1 <- data1[,subset.row,drop=FALSE]
        data2 <- data2[,subset.row,drop=FALSE]
    }

    scaling <- numeric(nrow(cell.vect))
    for (cell in seq_along(scaling)) {
        # For each cell, projecting both data sets onto the normalized correction vector for that cell.
        cur.cor.vect <- cell.vect[cell,]
        l2norm <- sqrt(sum(cur.cor.vect^2))
        cur.cor.vect <- cur.cor.vect/l2norm
        coords2 <- data2 %*% cur.cor.vect
        coords1 <- data1 %*% cur.cor.vect

        # Identifying the probability for that cell in its current data set, and the corresponding quantile in the other.
        prob2 <- (rank(coords2)[cell] - 0.5)/length(coords2)
        rank1 <- pmax(1, round(prob2 * length(coords1)))
        ord1 <- order(coords1)[rank1]
        
        # Adjusting the length of the correction vector so that the correction will match the quantiles.
        quan2 <- coords2[cell]
        quan1 <- coords1[ord1]
        scaling[cell] <- (quan1 - quan2)/l2norm
    }

    return(scaling * correction)
}

get.bio.span <- function(exprs, ndim, subset.row=NULL, pc.approx=FALSE)
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

    if (!pc.approx) { 
        S <- svd(centred, nu=ndim, nv=nv)
    } else {
        # Note: can't use centers=, as that subtracts from columns not rows.
        S <- irlba::irlba(centred, nu=ndim, nv=nv, work=ndim+20) 
    }
   
    # Returning the basis vectors directly if there was no subsetting. 
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
    S <- svd(t(A) %*% B)
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
            sintheta <- svd(t(A) - (t(A) %*% B) %*% t(B))$d[1]
        } else {
            sintheta <- svd(t(B) - (t(B) %*% A) %*% t(A))$d[1]
        }
        theta <- asin(min(1, sintheta)) 
    }
    
    list(angle=180*theta/pi, nshared=shared)
}

cosine.norm <- function(X, mode=c("matrix", "all", "l2norm"))
# Computes the cosine norm, with some protection from zero-length norms.
{
    mode <- match.arg(mode)
    out <- .Call(cxx_cosine_norm, X, mode!="l2norm")
    names(out) <- c("matrix", "l2norm")
    switch(mode, all=out, matrix=out$matrix, l2norm=out$l2norm)
}

bpl.get.knnx <- function(data, query, k, BPPARAM) 
# Splits up the query and searches for nearest neighbors in the data.
{
    nworkers <- bpworkers(BPPARAM)
    if (nworkers > 1) {
        by.core <- vector("list", nworkers)
        assignments <- cut(seq_len(nrow(query)), nworkers)
        for (i in seq_len(nworkers)) {
            by.core[[i]] <- query[assignments==levels(assignments)[i],,drop=FALSE]
        }
    } else {
        by.core <- list(query)
    }

    to.run <- bplapply(by.core, FUN=get.knnx, data=data, k=k, BPPARAM=BPPARAM)
    do.call(mapply, c(to.run, FUN=rbind, SIMPLIFY = FALSE))
}
