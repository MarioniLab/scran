#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame Rle
fastMNN <- function(..., k=20, cos.norm=TRUE, ndist=3, d=50, approximate=FALSE, 
    irlba.args=list(), subset.row=NULL, auto.order=FALSE, pc.input=FALSE,
    compute.variances=FALSE, assay.type="logcounts", get.spikes=FALSE, 
    BNPARAM=NULL, BPPARAM=SerialParam()) 
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

    # Creating the PCA input, if it is not already low-dimensional.
    if (!pc.input) {
        out <- .SCEs_to_matrices(batches, assay.type=assay.type, subset.row=subset.row, get.spikes=get.spikes)
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

    # Other pre-flight checks.
    var.kept <- rep(1, nbatches)
    re.order <- NULL 
    if (!is.logical(auto.order)) {
        re.order <- as.integer(auto.order)
        if (!identical(sort(re.order), seq_len(nbatches))) {
            stop("integer 'auto.order' must contain a permutation of 1:nbatches") 
        }
        auto.order <- FALSE
    }
    use.order <- !is.null(re.order)

    mnn.pairings <- vector("list", nbatches-1L)
    for (bdx in 2:nbatches) {

        if (auto.order) {
            # Automatically choosing to merge batches with the largest number of MNNs.
            if (bdx==2L) {
                d.out <- .define_first_merge(pc.mat, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
                refdata <- pc.mat[[d.out$first]]
                curdata <- pc.mat[[d.out$second]]
                mnn.sets <- d.out$pairs
                precomp <- d.out$precomputed
                processed <- c(d.out$first, d.out$second)
            } else {
                d.out <- .define_next_merge(refdata, pc.mat, processed, precomp, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
                curdata <- pc.mat[[d.out$other]]
                mnn.sets <- d.out$pairs
                processed <- c(processed, d.out$other)
            }
        } else {
            # Using the suggested merge order, or the supplied order in ...
            if (bdx==2L) {
                ref.idx <- if (use.order) re.order[1] else 1L
                refdata <- pc.mat[[ref.idx]]
                processed <- ref.idx
            }
            cur.idx <- if (use.order) re.order[bdx] else bdx
            curdata <- pc.mat[[cur.idx]]
            mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
            processed <- c(processed, cur.idx)
        }

        # Estimate the overall batch vector.
        ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        overall.batch <- colMeans(ave.out$averaged)

        # Remove variation along the batch vector.
        if (compute.variances) {
            var.before <- .compute_intra_var(refdata, curdata, pc.mat, processed)
        }

        refdata <- .center_along_batch_vector(refdata, overall.batch)
        curdata <- .center_along_batch_vector(curdata, overall.batch)

        if (compute.variances) {
            var.after <- .compute_intra_var(refdata, curdata, pc.mat, processed)
            var.kept[seq_len(bdx)] <- var.kept[seq_len(bdx)] * var.after/var.before
        }

        # Repeating the MNN discovery, now that the spread along the batch vector is removed.
#        re.mnn.sets <- find.mutual.nn(refdata, curdata, k1=k, k2=k, BPPARAM=BPPARAM)
#        re.ave.out <- .average_correction(refdata, re.mnn.sets$first, curdata, re.mnn.sets$second)

        # Recompute correction vectors and apply them.
        re.ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        curdata <- .tricube_weighted_correction(curdata, re.ave.out$averaged, re.ave.out$second, k=k, ndist=ndist, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

        mnn.pairings[[bdx-1L]] <- DataFrame(first=mnn.sets$first, second=mnn.sets$second + nrow(refdata))
        refdata <- rbind(refdata, curdata)
    }

    # Adjusting the output back to the specified order.
    if (auto.order || use.order) {
        ordering <- vector("list", nbatches)
        last <- 0L
        for (idx in processed) { 
            ncells <- nrow(pc.mat[[idx]])
            ordering[[idx]] <- last + seq_len(ncells)
            last <- last + ncells
        }
        ordering <- unlist(ordering)
        refdata <- refdata[ordering,,drop=FALSE]

        relocate <- ordering
        relocate[ordering] <- seq_along(ordering)
        for (x in seq_along(mnn.pairings)) {
            current <- mnn.pairings[[x]]
            current$first <- relocate[current$first]
            current$second <- relocate[current$second]
            mnn.pairings[[x]] <- current
        }
        
        var.kept[processed] <- var.kept
    }
    
    # Characterizing the original order of the batches and cells.
    if (!is.null(names(batches))) {
        batch.labels <- names(batches)
    } else {
        batch.labels <- seq_along(batches)
    }
    ncells <- vapply(pc.mat, FUN=nrow, FUN.VALUE=0L)
    batch.ids <- Rle(batch.labels, ncells)

    # Formatting the output.
    output <- list(corrected=refdata, batch=batch.ids, pairs=mnn.pairings, order=batch.labels[processed])
    if (!pc.input) {
        output$rotation <- metadata(pc.mat)$rotation
    }
    if (compute.variances) {
        output$lost.var <- 1 - var.kept
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
    corvec <- unname(corvec)/as.vector(npairs)
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

#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, BNPARAM=NULL, BPPARAM=SerialParam())
# Computing tricube-weighted correction vectors for individual cells,
# using the nearest neighbouring cells _involved in MNN pairs_.
{
    cur.uniq <- curdata[in.mnn,,drop=FALSE]
    safe.k <- min(k, nrow(cur.uniq))
    closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
    curdata + weighted.correction
}

############################################
# Variance calculation functions.

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colVars
.compute_intra_var <- function(reference, current, all.pcs, order) {
    all.var <- numeric(length(order))

    last <- 0L
    dm <- DelayedArray(reference)
    for (i in seq_len(length(order)-1L)) {
        cur.ncells <- nrow(all.pcs[[order[i]]])
        chosen <- last + seq_len(cur.ncells)
        all.var[i] <- sum(colVars(dm, rows=chosen))
        last <- last + cur.ncells
    }

    all.var[length(order)] <- sum(colVars(DelayedArray(current)))
    return(all.var)
}

############################################
# Auto-ordering functions.

#' @importFrom BiocNeighbors queryKNN buildNNIndex 
#' @importFrom BiocParallel SerialParam
.define_first_merge <- function(pc.mat, k=20, BNPARAM=NULL, BPPARAM=SerialParam())
# Find the pair of matrices in 'options' which has the greatest number of MNNs.
{
    precomputed <- lapply(pc.mat, buildNNIndex, BNPARAM=BNPARAM)
    max.pairs <- list()

    for (First in seq_along(precomputed)) {
        fdata <- pc.mat[[First]]
        for (Second in seq_len(First-1L)) {
            sdata <- pc.mat[[Second]]

            W21 <- queryKNN(BNINDEX=precomputed[[Second]], query=fdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
            W12 <- queryKNN(BNINDEX=precomputed[[First]], query=sdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
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

#' @importFrom BiocNeighbors queryKNN buildNNIndex 
#' @importFrom BiocParallel SerialParam
.define_next_merge <- function(refdata, pc.mat, processed, precomputed, k=20, BNPARAM=NULL, BPPARAM=SerialParam()) 
# Find the matrix in pc.mat[-processed] that has the greater number of MNNs to 'refdata'.
{
    max.pairs <- list()
    precomp.ref <- buildNNIndex(refdata, BNPARAM=BNPARAM)

    for (other in seq_along(pc.mat)[-processed]) {
        W21 <- queryKNN(BNINDEX=precomputed[[other]], query=refdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
        W12 <- queryKNN(BNINDEX=precomp.ref, query=pc.mat[[other]], k=k, BPPARAM=BPPARAM, get.distance=FALSE)
        out <- .Call(cxx_find_mutual_nns, W21$index, W12$index)
        names(out) <- c("first", "second")

        if (length(out$first) > length(max.pairs$first))  {
            max.pairs <- out
            max.other <- other
        }
    }

    return(list(other=max.other, pairs=max.pairs))
}

