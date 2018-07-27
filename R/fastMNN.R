#' @export
#' @importFrom BiocParallel SerialParam
fastMNN <- function(..., k=20, cos.norm=TRUE, ndist=3, d=50, approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, auto.order=FALSE, pc.input=FALSE,
    assay.type="logcounts", use.spikes=FALSE, BPPARAM=SerialParam()) 
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

    if (!pc.input) {
        out <- .SCEs_to_matrices(batches, assay.type=assay.type, subset.row=subset.row, use.spikes=use.spikes)
        batches <- out$batches
        subset.row <- out$subset.row

        if (!is.null(subset.row) && !identical(subset.row, seq_len(nrow(batches[[1]])))) { 
            batches <- lapply(batches, "[", i=subset.row, , drop=FALSE) # Need the extra comma!
        }
        if (cos.norm) { 
            batches <- lapply(batches, FUN=cosineNorm, mode="matrix")
        }

        pc.mat <- .multi_pca(batches, approximate=approximate, irlba.args=irlba.args, d=d, use.crossprod=TRUE)
    } else {
        .check_batch_consistency(batches, byrow=FALSE)
        pc.mat <- batches
    }

    mnn.pairings <- vector("list", nbatches-1L)
    for (bdx in 2:nbatches) {

        if (auto.order) {
            # Automatically choosing to merge batches with the largest number of MNNs.
            if (bdx==2L) {
                d.out <- .define_first_merge(pc.mat, k=k, BPPARAM=BPPARAM)
                refdata <- pc.mat[[d.out$first]]
                curdata <- pc.mat[[d.out$second]]
                mnn.sets <- d.out$pairs
                precomp <- d.out$precomputed
                processed <- c(d.out$first, d.out$second)
            } else {
                d.out <- .define_next_merge(refdata, pc.mat, processed, precomp, k=k, BPPARAM=BPPARAM)
                curdata <- pc.mat[[d.out$other]]
                mnn.sets <- d.out$pairs
                processed <- c(processed, d.out$other)
            }
        } else {
            # Otherwise, using consecutive batches in the supplied order.
            if (bdx==2L) {
                refdata <- pc.mat[[1]]
                processed <- 1L
            }
            curdata <- pc.mat[[bdx]]
            mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
            processed <- c(processed, bdx)
        }

        # Estimate the overall batch vector.
        ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        overall.batch <- colMeans(ave.out$averaged)

        # Remove variation along the batch vector.
        refdata <- .center_along_batch_vector(refdata, overall.batch)
        curdata <- .center_along_batch_vector(curdata, overall.batch)

        # Repeating the MNN discovery, now that the spread along the batch vector is removed.
#        re.mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
#        re.ave.out <- .average_correction(refdata, re.mnn.sets$first, curdata, re.mnn.sets$second)

        # Recompute correction vectors and apply them.
        re.ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        curdata <- .tricube_weighted_correction(curdata, re.ave.out$averaged, re.ave.out$second, k=k, ndist=ndist, BPPARAM=BPPARAM)

        mnn.pairings[[bdx-1L]] <- DataFrame(first=mnn.sets$first, second=mnn.sets$second + nrow(refdata))
        refdata <- rbind(refdata, curdata)
    }
    
    # Characterizing the original order of the batches and cells.
    if (!is.null(names(batches))) {
        batch.ids <- names(batches)
    } else {
        batch.ids <- seq_along(batches)
    }
    batch.ids <- batch.ids[processed]
    ncells <- vapply(pc.mat, FUN=nrow, FUN.VALUE=0L)[processed]

    batch.ids <- rep(batch.ids, ncells)
    cell.ids <- unlist(lapply(ncells, seq_len), use.names=FALSE)
    origin <- DataFrame(batch=batch.ids, cell=cell.ids)

    # Formatting the output.
    output <- list(corrected=refdata, origin=origin, pairs=mnn.pairings)
    if (!pc.input) {
        output$rotation <- metadata(pc.mat)$rotation
    }
    return(output)
}

############################################
# Correction-related functions.

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
# Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
# This removes any variation along the overall batch vector within each matrix.
{
    batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
    batch.loc <- as.vector(mat %*% batch.vec)
    central.loc <- mean(batch.loc)
    mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
    return(mat)
}

#' @importFrom kmknn queryKNN
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, BPPARAM=SerialParam())
# Computing tricube-weighted correction vectors for individual cells,
# using the nearest neighbouring cells _involved in MNN pairs_.
{
    cur.uniq <- curdata[in.mnn,,drop=FALSE]
    safe.k <- min(k, nrow(cur.uniq))
    closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BPPARAM=BPPARAM)

    middle <- ceiling(safe.k/2L)
    mid.dist <- closest$distance[,middle]
    rel.dist <- closest$distance / (mid.dist * ndist)
    rel.dist[rel.dist > 1] <- 1

    tricube <- (1 - rel.dist^3)^3
    weight <- tricube/rowSums(tricube)
    for (kdx in seq_len(safe.k)) {
        curdata <- curdata + correction[closest$index[,kdx],,drop=FALSE] * weight[,kdx]
    }
    
    return(curdata)
}

############################################
# Auto-ordering functions.

#' @importFrom kmknn queryKNN precluster
#' @importFrom BiocParallel SerialParam
.define_first_merge <- function(pc.mat, k=20, BPPARAM=SerialParam())
# Find the pair of matrices in 'options' which has the greatest number of MNNs.
{
    precomputed <- lapply(pc.mat, precluster)
    max.pairs <- list()

    for (First in seq_along(precomputed)) {
        fdata <- pc.mat[[First]]
        for (Second in seq_len(First-1L)) {
            sdata <- pc.mat[[Second]]

            W21 <- queryKNN(precomputed=precomputed[[Second]], query=fdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
            W12 <- queryKNN(precomputed=precomputed[[First]], query=sdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
            out <- .Call(cxx_find_mutual_nns, W21$index, W12$index)
            names(out) <- c("first", "second")

            if (length(out$first) > length(max.pairs$first))  {
                max.pairs <- out
                max.first <- First
                max.second <- Second
            }
        }
    }

    return(list(first=max.first, second=max.second, pairs=max.pairs, precomputed=precomputed))
}

#' @importFrom kmknn queryKNN precluster
#' @importFrom BiocParallel SerialParam
.define_next_merge <- function(refdata, pc.mat, processed, precomputed, k=20, BPPARAM=SerialParam()) 
# Find the matrix in pc.mat[-processed] that has the greater number of MNNs to 'refdata'.
{
    max.pairs <- list()
    precomp.ref <- precluster(refdata)

    for (other in seq_along(pc.mat)[-processed]) {
        W21 <- queryKNN(precomputed=precomputed[[other]], query=refdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
        W12 <- queryKNN(precomputed=precomp.ref, query=pc.mat[[other]], k=k, BPPARAM=BPPARAM, get.distance=FALSE)
        out <- .Call(cxx_find_mutual_nns, W21$index, W12$index)
        names(out) <- c("first", "second")

        if (length(out$first) > length(max.pairs$first))  {
            max.pairs <- out
            max.other <- other
        }
    }

    return(list(other=max.other, pairs=max.pairs))
}

