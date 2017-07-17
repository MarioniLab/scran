mnnCorrect <- function(..., inquiry.genes=NULL, hvg.genes=NULL, k=20, sigma=0.1, cos.norm=TRUE, svd.dim=2, order=NULL) 
# Performs correction based on the batches specified in the ellipsis.
#    
# written by Laleh Haghverdi
# with modifications by Aaron Lun
# created 7 April 2017
# last modified 22 June 2017
{
    batches <- batches0 <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { stop("at least two batches must be specified") }
   
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
    if (!is.null(inquiry.genes)) {
        batches0 <- lapply(batches0, "[", i=inquiry.genes, , drop=FALSE) # Need the extra comma!
    }
    if (!is.null(hvg.genes)) { 
        batches <- lapply(batches, "[", i=hvg.genes, , drop=FALSE)
    }
    inquiry.genes <- .subset_to_index(inquiry.genes, first, byrow=TRUE)
    hvg.genes <- .subset_to_index(hvg.genes, first, byrow=TRUE)
    inquiry.in.hvg <- inquiry.genes %in% hvg.genes 

    # Applying cosine normalization for MNN identification. 
    if (cos.norm) { batches <- lapply(batches, cosine.norm) }

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
    ref.batch <- batches[[ref]]
    ref.batch0 <- batches0[[ref]]
    
    num.mnn <- matrix(NA_integer_, nbatches, 2)
    output0 <- vector("list", nbatches)
    output0[[ref]] <- ref.batch0
    
    for (b in 2:nbatches) { 
        other.batch <- batches[[order[b]]]
        other.batch0 <- batches0[[order[b]]]
        
        # Finding pairs of mutual nearest neighbours.
        sets <- find.mutual.nn(ref.batch, other.batch, ref.batch0, other.batch0, k1=k, k2=k, sigma=sigma)
        s1 <- sets$set1
        s2 <- sets$set2

        if (svd.dim==0){
            correction <- t(sets$vect)
            correction0 <- t(sets$vect0)
        } else {
            ## Computing the biological subspace in both batches.
            ndim <- min(c(svd.dim, dim(ref.batch), dim(other.batch)))
            span1 <- get.bio.span(ref.batch0[,s1,drop=FALSE], inquiry.in.hvg, min(ndim, length(s1)))
            span2 <- get.bio.span(other.batch0[,s2,drop=FALSE], inquiry.in.hvg, min(ndim, length(s2)))
            #nshared <- find.shared.subspace(span1, span2, assume.orthonormal=TRUE, get.angle=FALSE)$nshared
            #if (nshared==0L) { warning("batches not sufficiently related") }
    
            #reduce the component in each span from the batch correction vector, span1 span2 order does not matter
            bv <- sets$vect
            bv0 <- sets$vect0      
            bio.comp <- bv0 %*% span1 %*% t(span1)
            correction0 <- t(bv0) - t(bio.comp)
            bio.comp <- t(correction0) %*% span2 %*% t(span2)
            correction0 <- correction0 - t(bio.comp)
            
            correction <- t(sets$vect)
        } 
        
        # Applying the correction and storing the numbers of nearest neighbors.
        other.batch <- other.batch + correction
        other.batch0 <- other.batch0 + correction0
        
        num.mnn[b,] <- c(length(s1), length(s2))
        output0[[b]] <- other.batch0

        # Expanding the reference batch to include the new, corrected data.
        ref.batch <- cbind(ref.batch, other.batch)
        ref.batch0 <- cbind(ref.batch0, other.batch0)
    }

    # Formatting output to be consistent with input.
    names(output0) <- names(batches0)
    list(corrected=output0, num.mnn=num.mnn)
}

find.mutual.nn <- function(exprs1, exprs2, exprs10, exprs20, k1, k2, sigma=0.1)
# Finds mutal neighbors between data1 and data2.
# Computes the batch correction vector for each cell in data2.
{
    data1 <- t(exprs1)
    data2 <- t(exprs2)
    
    data10 <- t(exprs10)
    data20 <- t(exprs20)
    
    n1 <- nrow(data1)
    n2 <- nrow(data2)
    n.total <- n1 + n2
   
    W21 <- FNN::get.knnx(data2, query=data1, k=k1)
    W12 <- FNN::get.knnx(data1, query=data2, k=k2)
    W <- sparseMatrix(i=c(rep(seq_len(n1), k1), rep(n1 + seq_len(n2), k2)),
                      j=c(n1 + W21$nn.index, W12$nn.index),
                      x=rep(1, n1*k1 + n2*k2), dims=c(n.total, n.total))

    W <- W * t(W) # elementwise multiplication to keep mutual nns only
    A <- which(W>0, arr.ind=TRUE) # row/col indices of mutual NNs
    set <- A

    # Computing the batch correction vector between MNN pairs.
    A1 <- A[,1]
    A1 <- A1[A1 <= n1]
    A2 <- A[,2] - n1
    A2 <- A2[A2 > 0]
    vect <- data1[A1,] - data2[A2,]    
    vect0 <- data10[A1,] - data20[A2,] 
    
    # Gaussian smoothing of individual correction vectors for MNN pairs.
    if (sigma==0) {
        G <- matrix(1, n2, n2)
    } else if (n2<3000) {
        dd2 <- as.matrix(dist(data2))
        G <- exp(-dd2^2/sigma)
    } else {
        kk <- min(length(A2),100)
        W <- get.knnx(data2[A2,], query=data2, k=kk)
        G <- matrix(0,n2,n2)
        for (i in seq_len(n2)) { 
            #G[i,A2[W$nn.index[i,]]]=W$}
            G[i,A2[W$nn.index[i,]]]=exp(-(W$nn.dist[i,])^2/sigma) 
        }
    }
    G <- (G+t(G)) /2
    
    #################
    D <- rowSums(G)
    nA2 <- tabulate(A2, nbins=n2)
    norm.dens <- t(G/(D*nA2))[,A2,drop=FALSE] # density normalized to avoid domination from dense parts
    batchvect <- norm.dens %*% vect 
    partitionf <- rowSums(norm.dens)
    partitionf[partitionf==0]<-1  # to avoid nans (instead get 0s)
    batchvect <- batchvect/partitionf

    batchvect0 <- norm.dens %*% vect0 
    batchvect0 <- batchvect0/partitionf
    
    # Report cells that are MNNs, and the correction vector per cell in data2.
    list(set1=unique(A1), set2=unique(A2), vect=batchvect, vect0=batchvect0)
}

get.bio.span <- function(exprs, inquiry.in.hvg, ndim) 
# Computes the basis matrix of the biological subspace of 'exprs'.
# The first 'ndim' dimensions are assumed to capture the biological subspace.
# Avoids extra dimensions dominated by technical noise, which will result in both 
# trivially large and small angles when using find.shared.subspace().
{
    keeph <- numeric(nrow(exprs))
    keeph[inquiry.in.hvg] <- 1
    exprs <- exprs * keeph
    exprs <- exprs - rowMeans(exprs) 
    S <- svd(exprs, nu=ndim, nv=0)
    S$u   
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

cosine.norm <- function(X)
# Computes the cosine norm, with some protection from zero-length norms.
{
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
}
