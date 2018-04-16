#' @importFrom scater librarySizeFactors normalize
#' @importFrom BiocGenerics sizeFactors
#' @importFrom S4Vectors split
#' @importFrom stats pt
.doublet_cluster <- function(x, clusters, size.factors=NULL, center.cluster.sf=TRUE, subset.row=NULL, threshold=0.05) 
# Finds evidence that a cluster is _not_ a doublet of two other clusters.
# Absence of such evidence should be a warning flag for that cluster.
# 
# written by Aaron Lun
# created 16 April 2018
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x)
    }

    # Centering the size factors.
    if (center.cluster.sf) {
        by.cluster <- split(seq_len(ncol(x)), clusters)
        for (i in by.cluster) {
            current <- size.factors[i]
            size.factors[i] <- current/mean(current)
        }
    }

    # Computing normalized counts.
    sce <- SingleCellExperiment(list(counts=x))
    sizeFactors(sce) <- size.factors
    sce <- normalize(sce, return_log=FALSE)

    # Computing mean, variance statistics off the normalized counts.
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    stats.out <- .get_var_stats(normcounts(sce), block=clusters, design=NULL, subset.row=subset.row)
    all.means <- stats.out$means
    all.vars <- stats.out$vars

    all.df <- stats.out$resid.df
    all.err <- sweep(all.vars, 2, all.df + 1, "/")
    all.df.comp <- sweep(all.err^2, 2, all.df, "/")
    
    clust.names <- colnames(all.means) 
    Nclusters <- length(clust.names)

    # Running through all pairs of clusters and testing against the third cluster.
    collected <- list()
    for (i1 in seq_len(Nclusters)) { 
        for (i2 in seq_len(i1-1L)) {
            
            others <- seq_len(Nclusters)[-c(i1, i2)]
            diffs <- all.means[,others,drop=FALSE] - (all.means[,i1] + all.means[,i2])
            errs <- all.err[,others,drop=FALSE] + all.err[,i1] + all.err[,i2]

            # Using Welch-Satterthwaite equation to get d.f. of combined error.
            df <- errs^2 / (all.df.comp[,others,drop=FALSE] + all.df.comp[,i1] + all.df.comp[,i2])

            # Computing t-statistics, p-value.
            t.stat <- diffs/sqrt(errs)
            pval <- pt(-abs(t.stat), df=df) * 2 

            # Computing BH adjusted p-values (i.e., combined Simes p-values).
            output.gene <- integer(length(others))
            output.min <- numeric(length(others))
            for (j in seq_along(others)) {
                adj.p <- p.adjust(pval[,j], method="BH")
                min.g <- which.min(adj.p)
                output.gene[j] <- min.g
                output.min[j] <- adj.p[min.g]
            }

            collected[[length(collected)+1]] <- DataFrame(
                current=clust.names[others],                        
                cluster1=clust.names[i1],
                cluster2=clust.names[i2],
                min.p=output.min,
                gene=subset.row[output.gene]
            )
        }
    }

    # Splitting the results by cluster.
    collected <- do.call(rbind, collected)
    host <- collected$current
    collected$current <- NULL
    output <- split(collected, host)

    # Returning results in order of largest p-value.
    for (i in seq_along(output)) {
        output[[i]] <- output[[i]][order(output[[i]]$min.p, decreasing=TRUE),]
    }
    return(output)
}
