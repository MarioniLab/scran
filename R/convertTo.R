# A function to convert between formats

setGeneric("convertTo", function(x, ...) standardGeneric("convertTo"))

setMethod("convertTo", "SCESet", function(x, type=c("edgeR", "DESeq2", "monocle"),
    fData.col=NULL, pData.col=NULL, ..., assay, use.all.sf=TRUE, normalize=TRUE, 
    subset.row=NULL, get.spikes=FALSE) {

    # Setting up the extraction.
    type <- match.arg(type)
    sf <- suppressWarnings(sizeFactors(x))
    if (type=="edgeR" || type=="DESeq2") { 
        fd <- fData(x)[,fData.col,drop=FALSE]
        pd <- pData(x)[,pData.col,drop=FALSE] 
    } else if (type=="monocle") {
        fd <- featureData(x)[,fData.col,drop=FALSE]
        pd <- phenoData(x)[,pData.col,drop=FALSE] 
    }
    if (missing(assay)) { assay <- "counts" }

    # Determining whether spikes should be retained.
    if (is.null(subset.row)) {
        subset.row <- .spike_subset(x, get.spikes)
    }
    if (!is.null(subset.row)) {
        subset.row <- .subset_to_index(subset.row, x)
    }

    # Collecting size factors for spikes.
    offset.index <- rep(1L, nrow(x))
    collected.sfs <- list(sf)
    if (is.null(sf)) { 
        collected.sfs[[1]] <- numeric(length(sf))
    } 
    for (st in whichSpike(x)) {
        spike.sf <- suppressWarnings(sizeFactors(x, type=st))
        if (!is.null(spike.sf)) {
            it <- length(collected.sfs) + 1L
            offset.index[isSpike(x, type=st)] <- it
            collected.sfs[[it]] <- spike.sf
        }
    }
    collected.offs <- do.call(rbind, collected.sfs)
    collected.offs <- log(collected.offs)
    collected.offs <- collected.offs - rowMeans(collected.offs)

    # Constructing objects of various types.
    if (type=="edgeR") {
        y <- DGEList(assayDataElement(x, assay), ...)
        if (ncol(fd)) { y$genes <- fd }
        if (!is.null(sf)) { 
            nf <- log(sf/y$samples$lib.size)
            nf <- exp(nf - mean(nf))
            y$samples$norm.factors <- nf
        }
        if (ncol(pd)) { y$samples <- cbind(y$samples, pd) }
        if (!is.null(subset.row)) { 
            y <- y[subset.row,] 
            offset.index <- offset.index[subset.row]            
        }
        if (!all(offset.index==1L) && use.all.sf) {
            y <- scaleOffset.DGEList(y, collected.offs[offset.index,,drop=FALSE])
        }
        return(y)

    } else if (type=="DESeq2") {
        dds <- DESeq2::DESeqDataSetFromMatrix(assayDataElement(x, assay), pd, ~1, ...)
        S4Vectors::mcols(dds) <- fd
        if (!is.null(sf)) { 
            sizeFactors(dds) <- sf
        }
        if (!is.null(subset.row)) { 
            dds <- dds[subset.row,] 
            offset.index <- offset.index[subset.row]            
        }
        if (!all(offset.index==1L) && use.all.sf) {
            DESeq2::normalizationFactors(dds) <- exp(collected.offs)[offset.index,,drop=FALSE]
        }
        return(dds)

    } else if (type=="monocle") {
        if (normalize) {
            if (is.null(sf)) { stop("size factors not defined for normalization") }
            cur.exprs <- t(t(counts(x))/sf)
            if (!all(offset.index==1L) && use.all.sf) {
                for (it in seq_along(collected.sfs)[-1]) {
                    current <- offset.index==it
                    cur.exprs[current,] <- t(t(counts(x)[current,,drop=FALSE])/collected.sfs[[it]])
                }
            }
        } else {
            cur.exprs <- assayDataElement(x, assay)
        }
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd, ...)
        if (!is.null(subset.row)) { out <- out[subset.row,] }
        return(out)

    }
})
