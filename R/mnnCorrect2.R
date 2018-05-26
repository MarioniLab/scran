#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom FNN get.knnx 
mnnCorrect2 <- function(..., k=20, sigma=0.1, cos.norm.in=TRUE, 
    d=50, approximate=FALSE, rand.seed=1000, irlba.args=list(), 
    subset.row=NULL, BPPARAM=SerialParam()) 
# Executes a faster and more accurate version of the MNN batch correction approach.
# 
# written by Aaron Lun
# created 26 May 2018
{

    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }
        
    prep.out <- prepare.input.data(batches, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.in, subset.row=subset.row)
    pc.mat <- .multi_pca(prep.out$In, d=d, approximate=approximate, rand.seed=rand.seed, irlba.args=irlba.args)

    refdata <- pc.mat[[1]]
    for (bdx in 2:nbatches) {
        curdata <- pc.mat[[bdx]]
        
        # Finding MNNs between batches and obtaining an estimate of the overall batch vector.
        mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
        ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        overall.batch <- colMeans(ave.out$averaged)

        # Projecting and _centering_ the length along the batch vector.
        refdata <- .center_along_batch_vector(refdata, overall.batch)
        curdata <- .center_along_batch_vector(curdata, overall.batch)

        # Repeating the MNN discovery, now that the spread along the batch vector is removed.
        re.mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
        re.ave.out <- .average_correction(refdata, re.mnn.sets$first, curdata, re.mnn.sets$second)
        corvec <- re.ave.out$averaged
        cur.uniq <- curdata[re.ave.out$second,,drop=FALSE]

        # Computing tricube-weighted correction vectors for individual cells.
        safe.k <- min(k, nrow(cur.uniq))
        closest <- get.knnx(query=curdata, data=cur.uniq, k=safe.k)
#        max.dist <- closest$nn.dist[,safe.k]
#        tricube <- (1 - (closest$nn.dist / max.dist)^3)^3
#        weight <- tricube/rowSums(tricube)
        weight <- 1/safe.k
        for (kdx in seq_len(safe.k)) {
            curdata <- curdata + corvec[closest$nn.index[,kdx],,drop=FALSE] * weight
        }

        # Combining this with the previous batch.
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
.multi_pca <- function(mat.list, d=50, approximate=FALSE, rand.seed=100, irlba.args=list()) 
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

    # rbinding all elements together for SVD.
    combined <- as.matrix(do.call(rbind, scaled))
    if (approximate) {
        set.seed(rand.seed)
        svd.out <- do.call(irlba::irlba, c(list(A=combined, nu=0, nv=d), irlba.args))
    } else {
        svd.out <- svd(combined, nu=0, nv=d)
    }

    # Projecting the scaled matrices back into this space.
    final <- centered
    for (idx in seq_along(centered)) {
        final[[idx]] <- crossprod(as.matrix(centered[[idx]]), svd.out$v)
    }
    return(final)
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
