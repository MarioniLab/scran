# A function to convert between formats

setGeneric("convertTo", function(x, ...) standardGeneric("convertTo"))

setMethod("convertTo", "SCESet", function(x, type=c("edgeR", "DESeq2", "monocle"),
    fData.col=NULL, pData.col=NULL, ..., assay, normalize=TRUE, 
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
        subset.row <- .spikeSubset(x, get.spikes)
    }
    if (!is.null(subset.row)) {
        subset.row <- .subset_to_index(subset.row, x)
    }

    if (type=="edgeR") {
        y <- DGEList(assayDataElement(x, assay), ...)
        if (ncol(fd)) { y$genes <- fd }
        if (!is.null(sf)) { 
            nf <- log(sf/y$samples$lib.size)
            nf <- exp(nf - mean(nf))
            y$samples$norm.factors <- nf
        }
        if (ncol(pd)) { y$samples <- cbind(y$samples, pd) }
        if (!is.null(subset.row)) { y <- y[subset.row,] }
        return(y)

    } else if (type=="DESeq2") {
        dds <- DESeq2::DESeqDataSetFromMatrix(assayDataElement(x, assay), pd, ~1, ...)
        S4Vectors::mcols(dds) <- fd
        sizeFactors(dds) <- sf
        if (!is.null(subset.row)) { dds <- dds[subset.row,] }
        return(dds)

    } else if (type=="monocle") {
        if (normalize) {
            if (is.null(sf)) { stop("size factors not defined for normalization") }
            cur.exprs <- t(t(counts(x))/sf)
        } else {
            cur.exprs <- assayDataElement(x, assay)
        }
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd, ...)
        if (!is.null(subset.row)) { out <- out[subset.row,] }
        return(out)

    }
})
