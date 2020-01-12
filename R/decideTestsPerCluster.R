#' Decide tests for each cluster
#'
#' Decide which tests (i.e., genes) are significant for differential expression between conditions in each cluster,
#' using the output of \code{\link{pseudoBulkDGE}}.
#' This mimics the \code{\link{decideTests}} functionality from \pkg{limma}.
#'
#' @param results A \linkS4class{List} containing the output of \code{\link{pseudoBulkDGE}}.
#' Each entry should be a DataFrame with the same number and order of rows,
#' containing at least a numeric \code{"PValue"} column (and usually a \code{"logFC"} column).
#' @param method String specifying whether the Benjamini-Hochberg correction should be applied across all clustesr
#' or separately within each cluster.
#' @param threshold Numeric scalar specifying the FDR threshold to consider genes as significant.
#' @param pval.field String containing the name of the column containing the p-value in each entry of \code{results}.
#' @param lfc.field String containing the name of the column containing the log-fold change.
#' Ignored if the column is not available Defaults to \code{"logFC"} if this field is available.
#'
#' @return
#' An integer matrix indicating whether each gene (row) is significantly DE between conditions for each cluster (column).
#' 
#' @details
#' If a log-fold change field is available and specified in \code{lfc.field}, values of \code{1}, \code{-1} and \code{0}
#' indicate that the gene is significantly upregulated, downregulated or not significant, respectively.
#' Note, the interpretation of \dQuote{up} and \dQuote{down} depends on the design and contrast in \code{\link{pseudoBulkDGE}}.
#'
#' Otherwise, if no log-fold change is available or if \code{lfc.field=NULL},
#' values of \code{1} or \code{0} indicate that a gene is significantly DE or not, respectively. 
#'
#' \code{NA} values indicate either that the relevant gene was low-abundance for a particular cluster and filtered out,
#' or that the DE comparison for that cluster was not possible (e.g., no residual d.f.).
#' 
#' @author Aaron Lun
#'
#' @examples
#' example(pseudoBulkDGE)
#' head(decideTestsPerCluster(out))
#' 
#' @seealso
#' \code{\link{pseudoBulkDGE}}, which generates the input to this function.
#'
#' \code{\link{decideTests}}, which inspired this function.
#'
#' @export
#' @importFrom stats p.adjust
decideTestsPerCluster <- function(results, method=c("separate", "global"), threshold=0.05, 
    pval.field="PValue", lfc.field="logFC") 
{
    method <- match.arg(method)
    all.p <- lapply(results, "[[", i="PValue")

    if (method=="separate") {
        all.p <- lapply(all.p, p.adjust, method="BH")
        all.p <- do.call(cbind, all.p)
        
    } else {
        all.p <- do.call(cbind, all.p)
        all.p[] <- p.adjust(all.p, method="BH")
    }        

    rownames(all.p) <- rownames(results[[1]])
    sig <- all.p <= threshold
    storage.mode(sig) <- "integer"

    if (!is.null(lfc.field) && !lfc.field %in% colnames(results[[1]])) {
        lfc.field <- NULL
    }
    if (!is.null(lfc.field)) {
        all.lfc <- do.call(cbind, lapply(results, "[[", i=lfc.field))
        sig <- sig * sign(all.lfc)
    }

    sig
}
