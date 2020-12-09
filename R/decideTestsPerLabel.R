#' Decide tests for each label
#'
#' Decide which tests (i.e., genes) are significant for differential expression between conditions in each label,
#' using the output of \code{\link{pseudoBulkDGE}}.
#' This mimics the \code{\link{decideTests}} functionality from \pkg{limma}.
#'
#' @param results A \linkS4class{List} containing the output of \code{\link{pseudoBulkDGE}}.
#' Each entry should be a DataFrame with the same number and order of rows,
#' containing at least a numeric \code{"PValue"} column (and usually a \code{"logFC"} column).
#'
#' For \code{summarizeTestsPerLabel}, this may also be a matrix produced by \code{decideTestsPerLabel}.
#' @param method String specifying whether the Benjamini-Hochberg correction should be applied across all clustesr
#' or separately within each label.
#' @param threshold Numeric scalar specifying the FDR threshold to consider genes as significant.
#' @param pval.field String containing the name of the column containing the p-value in each entry of \code{results}.
#' Defaults to \code{"PValue"}, \code{"P.Value"} or \code{"p.value"} based on fields in the first entry of \code{results}.
#' @param lfc.field String containing the name of the column containing the log-fold change.
#' Ignored if the column is not available Defaults to \code{"logFC"} if this field is available.
#' @param ... Further arguments to pass to \code{decideTestsPerLabel} if \code{results} is a List.
#'
#' @return
#' For \code{decideTestsPerLabel},
#' an integer matrix indicating whether each gene (row) is significantly DE between conditions for each label (column).
#'
#' For \code{summarizeTestsPerLabel},
#' an integer matrix containing the number of genes of each DE status (column) in each label (row).
#' 
#' @details
#' If a log-fold change field is available and specified in \code{lfc.field}, values of \code{1}, \code{-1} and \code{0}
#' indicate that the gene is significantly upregulated, downregulated or not significant, respectively.
#' Note, the interpretation of \dQuote{up} and \dQuote{down} depends on the design and contrast in \code{\link{pseudoBulkDGE}}.
#'
#' Otherwise, if no log-fold change is available or if \code{lfc.field=NULL},
#' values of \code{1} or \code{0} indicate that a gene is significantly DE or not, respectively. 
#'
#' \code{NA} values indicate either that the relevant gene was low-abundance for a particular label and filtered out,
#' or that the DE comparison for that label was not possible (e.g., no residual d.f.).
#' 
#' @author Aaron Lun
#'
#' @examples
#' example(pseudoBulkDGE)
#' head(decideTestsPerLabel(out))
#' summarizeTestsPerLabel(out)
#' 
#' @seealso
#' \code{\link{pseudoBulkDGE}}, which generates the input to this function.
#'
#' \code{\link{decideTests}}, which inspired this function.
#'
#' @export
#' @importFrom stats p.adjust
decideTestsPerLabel <- function(results, method=c("separate", "global"), threshold=0.05, 
    pval.field=NULL, lfc.field="logFC") 
{
    method <- match.arg(method)

    if (is.null(pval.field)) {
        pval.field <- intersect(c("PValue", "P.Value", "p.value"), colnames(results[[1]]))
        if (length(pval.field)==0) {
            stop("could not automatically determine 'pval.field'")
        }
        pval.field <- pval.field[1]
    }
    all.p <- lapply(results, "[[", i=pval.field)

    if (method=="separate") {
        all.p <- lapply(all.p, p.adjust, method="BH")
        all.p <- do.call(cbind, all.p)
        
    } else {
        all.p <- do.call(cbind, all.p)
        all.p[] <- p.adjust(all.p, method="BH")
    }        

    rownames(all.p) <- rownames(results[[1]])
    sig <- all.p <= threshold

    if (!is.null(lfc.field) && !lfc.field %in% colnames(results[[1]])) {
        lfc.field <- NULL
    }
    if (!is.null(lfc.field)) {
        all.lfc <- do.call(cbind, lapply(results, "[[", i=lfc.field))
        sig <- sig * sign(all.lfc)
    }

    storage.mode(sig) <- "integer"
    sig
}

#' @export
#' @rdname decideTestsPerLabel 
summarizeTestsPerLabel <- function(results, ...) {
    if (!is.matrix(results)) {
        results <- decideTestsPerLabel(results, ...)
    }

    output <- list()
    available <- sort(unique(as.vector(results)), na.last=TRUE)
    for (i in available) {
        output[[as.character(i)]] <- if (is.na(i)) {
            colSums(is.na(results)) 
        } else {
            colSums(results==i, na.rm=TRUE)
        }
    }

    do.call(cbind, output)
}
