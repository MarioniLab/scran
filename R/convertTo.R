#' Convert to other classes
#'
#' Convert a \linkS4class{SingleCellExperiment} object into other classes for entry into other analysis pipelines.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param type A string specifying the analysis for which the object should be prepared.
#' @param ... Other arguments to be passed to pipeline-specific constructors.
#' @param assay.type A string specifying which assay of \code{x} should be put in the returned object.
#' @param subset.row, get.spikes See \code{?"\link{scran-gene-selection}"}.
#' 
#' @return
#' For \code{type="edgeR"}, a DGEList object is returned containing the count matrix.
#' Size factors are converted to normalization factors.
#' Gene-specific \code{rowData} is stored in the \code{genes} element, and cell-specific \code{colData} is stored in the \code{samples} element.
#' 
#' For \code{type="DESeq2"}, a DESeqDataSet object is returned containing the count matrix and size factors.
#' Additional gene- and cell-specific data is stored in the \code{mcols} and \code{colData} respectively.
#' 
#' For \code{type="monocle"}, a CellDataSet object is returned containing the count matrix and size factors.
#' Additional gene- and cell-specific data is stored in the \code{featureData} and \code{phenoData} respectively.
#' 
#' @details
#' This function converts an SingleCellExperiment object into various other classes in preparation for entry into other analysis pipelines, as specified by \code{type}.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link[edgeR]{DGEList}},
#' \code{\link[DESeq2:DESeqDataSet]{DESeqDataSetFromMatrix}},
#' \code{\link[monocle]{newCellDataSet}}, for specific class constructors.
#' 
#' @examples
#' data(example.sce)
#' y <- example.sce
#' sizeFactors(y) <- 2^rnorm(ncol(y)) # Adding some additional embellishments.
#' rowData(y)$SYMBOL <- paste0("X", seq_len(nrow(y)))
#' y$other <- sample(LETTERS, ncells, replace=TRUE)
#' 
#' # Converting to various objects.
#' convertTo(y, type="edgeR")
#' convertTo(y, type="DESeq2")
#' convertTo(y, type="monocle")
#' 
#' @export
#' @importFrom BiocGenerics sizeFactors as.data.frame "sizeFactors<-"
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom SingleCellExperiment isSpike
#' @importFrom S4Vectors "mcols<-" 
#' @importFrom edgeR DGEList "[.DGEList" scaleOffset.DGEList
convertTo <- function(x, type=c("edgeR", "DESeq2", "monocle"), ..., assay.type=1, subset.row=NULL) {
    type <- match.arg(type)
    if (type=="edgeR" || type=="DESeq2") { 
        fd <- rowData(x)
        pd <- colData(x)
    } else if (type=="monocle") {
        fd <- Biobase::AnnotatedDataFrame(as.data.frame(rowData(x), row.names=rownames(x)))
        pd <- Biobase::AnnotatedDataFrame(as.data.frame(colData(x)))
    }

    sf <- suppressWarnings(sizeFactors(x))
    subset.row <- .subset_to_index(subset.row, x)

    # Constructing objects of various types.
    if (type=="edgeR") {
        y <- DGEList(assay(x, i=assay.type)[subset.row,,drop=FALSE], ...)
        if (ncol(fd)) { 
            y$genes <- fd[subset.row,,drop=FALSE]
        }

        if (!is.null(sf)) { 
            nf <- log(sf/y$samples$lib.size)
            nf <- exp(nf - mean(nf))
            y$samples$norm.factors <- nf
        }

        if (ncol(pd)) { 
            y$samples <- cbind(y$samples, pd) 
        }
        return(y)

    } else if (type=="DESeq2") {
        dds <- DESeq2::DESeqDataSetFromMatrix(assay(x, i=assay.type)[subset.row,,drop=FALSE], pd, ~1, ...)
        rowData(dds) <- fd[subset.row,]
        if (!is.null(sf)) { 
            sizeFactors(dds) <- sf
        }
        return(dds)

    } else if (type=="monocle") {
        cur.exprs <- assay(x, i=assay.type)[subset.row,,drop=FALSE]
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd[subset.row,], ...)
        if (!is.null(sf)) { 
            sizeFactors(out) <- sf 
        }
        return(out)
    }
}
