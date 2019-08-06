#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
correlateGenes <- function(stats) 
# Combine statistics for gene pair correlations into a single statistic per gene.
#
# written by Aaron Lun
# created 5 January 2019
{
    pool <- union(stats$gene1, stats$gene2)
    m1 <- match(stats$gene1, pool)
    m2 <- match(stats$gene2, pool)
    by.gene <- combine_rho(length(pool), m1 - 1L, m2 - 1L,
        stats$rho, stats$p.value, stats$limited, order(stats$p.value) - 1L)
    
    out <- DataFrame(gene=pool, rho=by.gene[[2]], p.value=by.gene[[1]], 
        FDR=p.adjust(by.gene[[1]], method="BH"), limited=by.gene[[3]])
    .is_sig_limited(out)

    out
}
