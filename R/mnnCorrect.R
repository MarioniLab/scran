mnnCorrect <- function(..., k=20, sigma=0.1, cos.norm=TRUE, svd.dim=20, subset.row=NULL, order=NULL, 
                       pc.approx=FALSE, exact.kernel=TRUE, kernel.k=100, BPPARAM=SerialParam())
# Performs correction based on the batches specified in the ellipsis.
#    
# written by Laleh Haghverdi
# with modifications by Aaron Lun
# created 7 April 2017
# last modified 3 October 2017
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
    if (!is.null(subset.row)) { 
        batches <- lapply(batches, "[", i=subset.row, , drop=FALSE) # Need the extra comma!
    }

    # Applying cosine normalization for MNN identification. 
    if (cos.norm) { 
        batches <- lapply(batches, cosine.norm) 
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
    ref.batch <- t(batches[[ref]])

    num.mnn <- matrix(NA_integer_, nbatches, 2)
    output <- vector("list", nbatches)
    output[[ref]] <- batches[[ref]]
    
    for (b in 2:nbatches) { 
        other.batch <- t(batches[[order[b]]])
        
        # Finding pairs of mutual nearest neighbours.
        sets <- find.mutual.nn(ref.batch, other.batch, k1=k, k2=k, BPPARAM=BPPARAM)
        s1 <- sets$first
        s2 <- sets$second      

        kernel <- construct.smoothing.kernel(other.batch, sigma=sigma, exact=exact.kernel,
                                             kk=kernel.k, mnn.set=s2, BPPARAM=BPPARAM)
        correction <- compute.correction.vectors(ref.batch, other.batch, s1, s2, kernel)

        if (!is.na(svd.dim)){
            # Computing the biological subspace in both batches.
            ndim <- min(c(svd.dim, dim(ref.batch), dim(other.batch)))
            span1 <- get.bio.span(t(ref.batch[s1,,drop=FALSE]), ndim=min(ndim, length(s1)), pc.approx=pc.approx)
            span2 <- get.bio.span(t(other.batch[s2,,drop=FALSE]), ndim=min(ndim, length(s2)), pc.approx=pc.approx)

            #nshared <- find.shared.subspace(span1, span2, assume.orthonormal=TRUE, get.angle=FALSE)$nshared
            #if (nshared==0L) { warning("batches not sufficiently related") }
    
            # Reduce the component in each span from the batch correction vector.
            correction <- subtract.bio(correction, span1, span2)
        } 
        
        # Applying the correction and storing the numbers of nearest neighbors.
        other.batch <- other.batch + correction       
        num.mnn[b,] <- c(length(s1), length(s2))
        output[[b]] <- t(other.batch)

        # Expanding the reference batch to include the new, corrected data.
        ref.batch <- rbind(ref.batch, other.batch)
    }

    # Formatting output to be consistent with input.
    names(output) <- names(batches)
    list(corrected=output, num.mnn=num.mnn)
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

construct.smoothing.kernel <- function(data, sigma=0.1, exact=TRUE, kk=100, mnn.set, BPPARAM) 
# Constructs a Gaussian smoothing kernel, using all distances or the closest 100 cells.
{ 
    N <- nrow(data)
    if (is.na(sigma)) {
        G <- matrix(1, N, N)
    } else if (exact) { 
        dd2 <- as.matrix(dist(data))
        G <- exp(-dd2^2/sigma)
    } else {
        mnn.set <- unique(mnn.set)
        kk <- min(length(mnn.set), kk)
        W <- bpl.get.knnx(data[mnn.set,,drop=FALSE], query=data, k=kk, BPPARAM=BPPARAM)
        vals <- as.vector(exp(-W$nn.dist^2/sigma))
        i.dex <- rep(seq_len(N), kk)
        j.dex <- as.vector(mnn.set[W$nn.index])
        G <- sparseMatrix(i=i.dex, j=j.dex, x=vals, dims=c(N, N), dimnames=NULL)
    }
    return(G)
}

compute.correction.vectors <- function(data1, data2, mnn1, mnn2, kernel) 
# Computes the batch correction vector for each cell in data2.
{
    vect <- data1[mnn1,] - data2[mnn2,]    
       
    # Density normalized to avoid domination from dense parts
    D <- rowSums(kernel)
    nA2 <- tabulate(mnn2, nbins=nrow(data2))
    norm.dens <- t(kernel/(D*nA2))[,mnn2,drop=FALSE]

    # Computing normalized batch correction vectors.
    batchvect <- norm.dens %*% vect 
    partitionf <- rowSums(norm.dens)
    partitionf[partitionf==0] <- 1  # to avoid nans (instead get 0s)
    return(batchvect/partitionf)
}

get.bio.span <- function(exprs, ndim, pc.approx=FALSE)
# Computes the basis matrix of the biological subspace of 'exprs'.
# The first 'ndim' dimensions are assumed to capture the biological subspace.
# Avoids extra dimensions dominated by technical noise, which will result in both 
# trivially large and small angles when using find.shared.subspace().
{
    centred <- exprs - rowMeans(exprs)
    centred <- as.matrix(centred)
    if (!pc.approx) { 
        S <- svd(centred, nu=ndim, nv=0)
    } else {
        S <- irlba::irlba(centred, nu=ndim, nv=0, work=ndim+20) # can't use centers=, as that subtracts from columns not rows.
    }
    return(S$u)
}

subtract.bio <- function(correction, span1, span2) 
# Subtracts the biological basis vectors from the correction vectors.
# Note that the order of span1 and span2 does not matter.
{ 
    bio.comp <- correction %*% span1 %*% t(span1)
    correction <- correction - bio.comp
    bio.comp <- correction %*% span2 %*% t(span2)
    correction <- correction - bio.comp
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

cosine.norm <- function(X)
# Computes the cosine norm, with some protection from zero-length norms.
{
    .Call(cxx_cosine_norm, X)
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
