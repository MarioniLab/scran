#' Identify HVGs
#' 
#' Define a set of highly variable genes, based on variance modelling statistics
#' from \code{\link{modelGeneVar}} or related functions.
#'
#' @param stats A \linkS4class{DataFrame} of variance modelling statistics with one row per gene.
#' Alternatively, a \linkS4class{SummarizedExperiment} object, in which case it is supplied to \code{\link{modelGeneVar}} to generate the required DataFrame.
#' @param var.field String specifying the column of \code{stats} containing the relevant metric of variation.
#' @param n Integer scalar specifying the number of top HVGs to report.
#' @param prop Numeric scalar specifying the proportion of genes to report as HVGs.
#' @param var.threshold Numeric scalar specifying the minimum threshold on the metric of variation.
#' @param fdr.field String specifying the column of \code{stats} containing the adjusted p-values.
#' If \code{NULL}, no filtering is performed on the FDR.
#' @param fdr.threshold Numeric scalar specifying the FDR threshold.
#' @param row.names Logical scalar indicating whether row names should be reported.
#'
#' @return 
#' A character vector containing the names of the most variable genes, if \code{row.names=TRUE}.
#'
#' Otherwise, an integer vector specifying the indices of \code{stats} containing the most variable genes.
#'
#' @details
#' This function will identify all genes where the relevant metric of variation is greater than \code{var.threshold}.
#' By default, this means that we retain all genes with positive values in the \code{var.field} column of \code{stats}.
#' If \code{var.threshold=NULL}, the minimum threshold on the value of the metric is not applied.
#' 
#' If \code{fdr.threshold} is specified, we further subset to genes that have FDR less than or equal to \code{fdr.threshold}.
#' By default, FDR thresholding is turned off as \code{\link{modelGeneVar}} and related functions 
#' determine significance of large variances \emph{relative} to other genes.
#' This can be overly conservative if many genes are highly variable.
#'
#' If \code{n=NULL} and \code{prop=NULL}, the resulting subset of genes is directly returned.
#' Otherwise, the top set of genes with the largest values of the variance metric are returned,
#' where the size of the set is defined as the larger of \code{n} and \code{prop*nrow(stats)}.
#' 
#' @seealso
#' \code{\link{modelGeneVar}} and friends, to generate \code{stats}.
#'
#' \code{\link{modelGeneCV2}} and friends, to also generate \code{stats}.
#' 
#' @author Aaron Lun
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' stats <- modelGeneVar(sce)
#' str(getTopHVGs(stats))
#' str(getTopHVGs(stats, fdr.threshold=0.05)) # more stringent
#'
#' # Or directly pass in the SingleCellExperiment:
#' str(getTopHVGs(sce))
#'
#' # Alternatively, use with the coefficient of variation:
#' stats2 <- modelGeneCV2(sce)
#' str(getTopHVGs(stats2, var.field="ratio"))
#' 
#' @export
#' @importFrom utils head
getTopHVGs <- function(stats, var.field="bio", n=NULL, prop=NULL, var.threshold=0,
    fdr.field="FDR", fdr.threshold=NULL, row.names=!is.null(rownames(stats))) 
{
    if (is(stats, "SummarizedExperiment")) {
        stats <- modelGeneVar(stats)
    }

    survivors <- seq_len(nrow(stats))

    if (!is.null(fdr.threshold)) {
        fdr <- stats[[fdr.field]]
        keep <- !is.na(fdr) & fdr <= fdr.threshold
        survivors <- survivors[keep]
        stats <- stats[keep,,drop=FALSE]
    }

    if (!is.null(var.threshold)) {
        var <- stats[[var.field]]
        keep <- !is.na(var) & var > var.threshold
        survivors <- survivors[keep]
        stats <- stats[keep,,drop=FALSE]
    }

    o <- order(stats[[var.field]], decreasing=TRUE)
    if (!is.null(n) || !is.null(prop)) {
        n <- max(n, round(prop*nrow(stats)))
        o <- head(o, n)
    }

    if (row.names) {
        rownames(stats)[o]
    } else {
        survivors[o]
    }
}
