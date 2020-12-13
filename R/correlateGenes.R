#' Per-gene correlation statistics
#' 
#' Compute per-gene correlation statistics by combining results from gene pair correlations.
#' 
#' @param stats A \linkS4class{DataFrame} of pairwise correlation statistics, returned by \code{\link{correlatePairs}}.
#' 
#' @return
#' A \linkS4class{DataFrame} with one row per unique gene in \code{stats} and containing the fields:
#' \describe{
#' \item{\code{gene}:}{A field of the same type as \code{stats$gene1} specifying the gene identity.}
#' \item{\code{rho}:}{Numeric, the correlation with the largest magnitude across all gene pairs involving the corresponding gene.}
#' \item{\code{p.value}:}{Numeric, the Simes p-value for this gene.}
#' \item{\code{FDR}:}{Numeric, the adjusted \code{p.value} across all rows.}
#' }
#' 
#' @details
#' For each gene, all of its pairs are identified and the corresponding p-values are combined using Simes' method.
#' This tests whether the gene is involved in significant correlations to \emph{any} other gene.
#' Per-gene statistics are useful for identifying correlated genes without regard to what they are correlated with (e.g., during feature selection).
#' 
#' @seealso
#' \code{\link{correlatePairs}}, to compute \code{stats}.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' pairs <- correlatePairs(sce, iters=1e5, subset.row=1:100)
#'
#' g.out <- correlateGenes(pairs)
#' head(g.out)
#' 
#' @references
#' Simes RJ (1986).
#' An improved Bonferroni procedure for multiple tests of significance.
#' \emph{Biometrika} 73:751-754.
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
correlateGenes <- function(stats) {
    pool <- union(stats$gene1, stats$gene2)
    m1 <- match(stats$gene1, pool)
    m2 <- match(stats$gene2, pool)
    by.gene <- combine_rho(length(pool), m1 - 1L, m2 - 1L, stats$rho, stats$p.value, order(stats$p.value) - 1L)
    DataFrame(gene=pool, rho=by.gene[[2]], p.value=by.gene[[1]], FDR=p.adjust(by.gene[[1]], method="BH"))
}
