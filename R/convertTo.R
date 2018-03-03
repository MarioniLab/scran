#' @importFrom BiocGenerics sizeFactors as.data.frame "sizeFactors<-"
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom SingleCellExperiment isSpike
#' @importFrom S4Vectors "mcols<-" 
#' @importFrom edgeR DGEList "[.DGEList" scaleOffset.DGEList
.convert_to <- function(x, type=c("edgeR", "DESeq2", "monocle"),
                        row.fields=NULL, col.fields=NULL, ..., assay.type, use.all.sf=TRUE, normalize=TRUE, 
                        subset.row=NULL, get.spikes=FALSE) 
# A function to convert between formats
#
# written by Aaron Lun
# a while ago    
{
    # Setting up the extraction.
    type <- match.arg(type)
    sf <- suppressWarnings(sizeFactors(x))
    if (type=="edgeR" || type=="DESeq2") { 
        fd <- rowData(x)[,row.fields,drop=FALSE]
        pd <- colData(x)[,col.fields,drop=FALSE] 
    } else if (type=="monocle") {
        fd <- Biobase::AnnotatedDataFrame(as.data.frame(rowData(x)[,row.fields,drop=FALSE], row.names=rownames(x)))
        pd <- Biobase::AnnotatedDataFrame(as.data.frame(colData(x)[,col.fields,drop=FALSE]))
    }
    if (missing(assay.type)) { assay.type <- "counts" }

    # Determining which genes should be extracted.
    subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)

    # Collecting size factors for spikes.
    offset.index <- rep(1L, nrow(x))
    collected.sfs <- list(sf)
    if (is.null(sf)) { 
        collected.sfs[[1]] <- numeric(length(sf))
    } 
    for (st in spikeNames(x)) {
        spike.sf <- sizeFactors(x, type=st)
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
        y <- DGEList(assay(x, i=assay.type), ...)
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
        dds <- DESeq2::DESeqDataSetFromMatrix(assay(x, i=assay.type), pd, ~1, ...)
        mcols(dds) <- fd
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
            dimnames(cur.exprs) <- dimnames(x)
        } else {
            cur.exprs <- assay(x, i=assay.type)
        }
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd, ...)
        if (!is.null(subset.row)) { out <- out[subset.row,] }
        return(out)

    }
}

#' @export
setGeneric("convertTo", function(x, ...) standardGeneric("convertTo"))

#' @export
setMethod("convertTo", "SingleCellExperiment", .convert_to)

