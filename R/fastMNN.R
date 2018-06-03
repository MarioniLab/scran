#' @export
#' @importFrom BiocParallel SerialParam
fastMNN <- function(..., k=20, cos.norm=TRUE, d=50, ndist=3, approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, BPPARAM=SerialParam()) 
# A faster version of the MNN batch correction approach.
# 
# written by Aaron Lun
# created 26 May 2018
{

    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }
        
    prep.out <- prepare.input.data(batches, cos.norm.in=cos.norm, cos.norm.out=cos.norm, subset.row=subset.row)
    pc.mat <- .multi_pca(prep.out$In, d=d, approximate=approximate, irlba.args=irlba.args)

    refdata <- pc.mat[[1]]
    for (bdx in 2:nbatches) {
        curdata <- pc.mat[[bdx]]
        
        # Finding MNNs between batches and obtaining an estimate of the overall batch vector.
        mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
        ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        overall.batch <- colMeans(ave.out$averaged)

        # Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
        refdata <- .center_along_batch_vector(refdata, overall.batch)
        curdata <- .center_along_batch_vector(curdata, overall.batch)

        # Repeating the MNN discovery, now that the spread along the batch vector is removed.
#        re.mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
#        re.ave.out <- .average_correction(refdata, re.mnn.sets$first, curdata, re.mnn.sets$second)

        # Recomputing the correction vectors after removing the within-batch variation.
        re.ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        
        curdata <- .tricube_weighted_correction(curdata, re.ave.out$averaged, re.ave.out$second, k=k, ndist=ndist)
        refdata <- rbind(refdata, curdata)
    }
    
    # Figuring out what output to return.
    ncells <- vapply(batches, FUN=ncol, FUN.VALUE=0L)
    if (!is.null(names(batches))) {
        batch.ids <- rep(names(batches), ncells)
    } else {
        batch.ids <- rep(seq_along(ncells), ncells)
    }
    return(list(corrected=refdata, batch=batch.ids))
}
 
#' @importFrom DelayedArray DelayedArray t
#' @importFrom DelayedMatrixStats rowMeans2
.multi_pca <- function(mat.list, d=50, approximate=FALSE, irlba.args=list()) 
# This function performs a multi-sample PCA, weighting the contribution of each 
# sample to the gene-gene covariance matrix to avoid domination by samples with
# a large number of cells. Expects cosine-normalized and subsetted expression matrices.
{
    all.centers <- 0
    for (idx in seq_along(mat.list)) {
        current <- DelayedArray(mat.list[[idx]])
        centers <- rowMeans2(current)
        all.centers <- all.centers + centers
        mat.list[[idx]] <- current
    }
    all.centers <- all.centers/length(mat.list) # grand average of centers (not using batch-specific centers, which makes compositional assumptions).

    centered <- scaled <- mat.list
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        current <- current - centers # centering each batch by the grand average.
        centered[[idx]] <- current
        current <- current/sqrt(ncol(current)) # downweighting samples with many cells.
        scaled[[idx]] <- t(current)
    }

    # Performing an SVD.
    if (approximate) {
        svd.out <- .fast_irlba(scaled, nv=d, irlba.args=irlba.args)
    } else {
        combined <- as.matrix(do.call(rbind, scaled))
        svd.out <- svd(combined, nu=0, nv=d)
    }

    # Projecting the scaled matrices back into this space.
    final <- centered
    for (idx in seq_along(centered)) {
        final[[idx]] <- crossprod(as.matrix(centered[[idx]]), svd.out$v)
    }
    return(final)
}

#' @importFrom DelayedArray t
.fast_irlba <- function(mat.list, nv, irlba.args=list()) 
# Performs a quick irlba by performing the SVD on XtX or XXt,
# and then obtaining the V vector from one or the other.
{
    nrows <- sum(vapply(mat.list, FUN=nrow, FUN.VALUE=0L))
    ncols <- ncol(mat.list[[1]])
    
    # Creating the cross-product without actually rbinding the matrices.
    # This avoids creating a large temporary matrix.
    flipped <- nrows > ncols
    if (flipped) {
        final <- matrix(0, ncols, ncols)
        for (mdx in seq_along(mat.list)) {
            curmat <- as.matrix(mat.list[[mdx]])
            final <- final + crossprod(curmat)
        }

    } else {
        final <- matrix(0, nrows, nrows)
        last1 <- 0L
        for (first in seq_along(mat.list)) {
            fmat <- as.matrix(mat.list[[first]])
            fdx <- last1 + seq_len(nrow(fmat))
            
            last2 <- 0L
            for (second in seq_along(mat.list)) {
                smat <- mat.list[[second]]
                sdx <- last2 + seq_len(nrow(smat))
                final[fdx,sdx] <- as.matrix(fmat %*% t(smat))
                last2 <- last2 + nrow(smat)
            }
            
            last1 <- last1 + nrow(fmat)
        }
    }

    svd.out <- do.call(irlba::irlba, c(list(A=final, nv=nv, nu=0), irlba.args))
    svd.out$d <- sqrt(svd.out$d)

    if (flipped) {
        # XtX means that the V in the irlba output is the original Vm,
        # which can be directly returned.
        return(svd.out)
    }

    # Otherwise, XXt means that the V in the irlba output is the original U.
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

.average_correction <- function(refdata, mnn1, curdata, mnn2)
# Computes correction vectors for each MNN pair, and then
# averages them for each MNN-involved cell in the second batch.
{
    corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
    corvec <- rowsum(corvec, mnn2)
    npairs <- table(mnn2)
    stopifnot(identical(names(npairs), rownames(corvec)))
    corvec <- corvec/as.vector(npairs)
    list(averaged=corvec, second=as.integer(names(npairs)))
}

.center_along_batch_vector <- function(mat, batch.vec) 
# This removes any variation along the overall batch vector within each matrix.
{
    batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
    batch.loc <- as.vector(mat %*% batch.vec)
    central.loc <- mean(batch.loc)
    mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
    return(mat)
}

#' @importFrom FNN get.knnx 
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3)
# Computing tricube-weighted correction vectors for individual cells,
# using the nearest neighbouring cells _involved in MNN pairs_.
{
    cur.uniq <- curdata[in.mnn,,drop=FALSE]
    safe.k <- min(k, nrow(cur.uniq))
    closest <- get.knnx(query=curdata, data=cur.uniq, k=safe.k)

    middle <- ceiling(safe.k/2L)
    mid.dist <- closest$nn.dist[,middle]
    rel.dist <- closest$nn.dist / (mid.dist * ndist)
    rel.dist[rel.dist > 1] <- 1

    tricube <- (1 - rel.dist^3)^3
    weight <- tricube/rowSums(tricube)
    for (kdx in seq_len(safe.k)) {
        curdata <- curdata + correction[closest$nn.index[,kdx],,drop=FALSE] * weight[,kdx]
    }
    
    return(curdata)
}
