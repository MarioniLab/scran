# A function to convert between formats

setGeneric("convertTo", function(x, ...) standardGeneric("convertTo"))

setMethod("convertTo", "SCESet", function(x, type=c("edgeR", "DESeq2", "monocle"),
    fData.col=NULL, pData.col=NULL, ..., assay, normalize=TRUE, get.spikes=FALSE) {

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

    if (type=="edgeR") {
        y <- DGEList(assayDataElement(x, assay), ...)
        if (ncol(fd)) { y$genes <- fd }
        if (!is.null(sf)) { 
            nf <- log(sf/y$samples$lib.size)
            nf <- exp(nf - mean(nf))
            y$samples$norm.factors <- nf
        }
        if (ncol(pd)) { y$samples <- cbind(y$samples, pd) }
        if (!get.spikes) { y <- y[!isSpike(x),] }
        return(y)

    } else if (type=="DESeq2") {
        dds <- DESeq2::DESeqDataSetFromMatrix(assayDataElement(x, assay), pd, ~1, ...)
        S4Vectors::mcols(dds) <- fd
        sizeFactors(dds) <- sf
        if (!get.spikes) { dds <- dds[!isSpike(x),] }
        return(dds)

    } else if (type=="monocle") {
        cur.exprs <- assayDataElement(x, assay)
        if (normalize) {
            if (is.null(sf)) { stop("size factors not defined for normalization") }
            cur.exprs <- cpm.default(cur.exprs, lib.size=sizeFactors(x)*1e6, log=FALSE, prior.count=0)
        }
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd, ...)
        if (!get.spikes) { out <- out[!isSpike(x),] }
        return(out)

    }
})
