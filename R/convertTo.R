#' Convert to other classes
#'
#' Convert a \linkS4class{SingleCellExperiment} object into other classes for entry into other analysis pipelines.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param type A string specifying the analysis for which the object should be prepared.
#' @param ... Other arguments to be passed to pipeline-specific constructors.
#' @param assay.type A string specifying which assay of \code{x} should be put in the returned object.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
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
#' library(scuttle)
#' sce <- mockSCE()
#' 
#' # Adding some additional embellishments.
#' sizeFactors(sce) <- 2^rnorm(ncol(sce)) 
#' rowData(sce)$SYMBOL <- paste0("X", seq_len(nrow(sce)))
#' sce$other <- sample(LETTERS, ncol(sce), replace=TRUE)
#' 
#' # Converting to various objects.
#' convertTo(sce, type="edgeR")
#' convertTo(sce, type="DESeq2")
#' 
#' @export
#' @importFrom BiocGenerics sizeFactors as.data.frame "sizeFactors<-"
#' @importFrom SummarizedExperiment rowData colData assay rowData<-
#' @importFrom S4Vectors "mcols<-" 
#' @importFrom edgeR DGEList "[.DGEList" scaleOffset.DGEList
#' @importFrom scuttle .subset2index
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
    subset.row <- .subset2index(subset.row, x)

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
        rowData(dds) <- fd[subset.row,,drop=FALSE]
        if (!is.null(sf)) { 
            sizeFactors(dds) <- sf
        }
        return(dds)

    } else if (type=="monocle") {
        cur.exprs <- assay(x, i=assay.type)[subset.row,,drop=FALSE]
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd[subset.row,,drop=FALSE], ...)
        if (!is.null(sf)) { 
            sizeFactors(out) <- sf 
        }
        return(out)
    }
}
