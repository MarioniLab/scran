#' @importFrom scater librarySizeFactors normalize
#' @importFrom BiocGenerics "sizeFactors<-"
#' @importFrom stats p.adjust
.doublet_cluster <- function(x, clusters, subset.row=NULL, threshold=0.05) 
# Finds evidence that a cluster is _not_ a doublet of two other clusters.
# Absence of such evidence should be a warning flag for that cluster.
# 
# written by Aaron Lun
# created 16 April 2018
{
    # Computing normalized counts using the library size (looking for compositional differences!)
    sce <- SingleCellExperiment(list(counts=x))
    sizeFactors(sce) <- librarySizeFactors(x)
    sce <- normalize(sce, return_log=TRUE)
    degs <- findMarkers(sce, clusters=clusters, subset.row=subset.row, full.stats=TRUE)

    # Running through all pairs of clusters and testing against the third cluster.
    collected <- list()
    all.clusters <- names(degs)
    for (ref in all.clusters) {
        ref.stats <- degs[[ref]]
        remnants <- setdiff(all.clusters, ref)
        best.N <- Inf
        best.p <- 1
        best.gene <- NULL
        best.parents <- NULL

        for (i1 in seq_along(remnants)) {
            stats1 <- ref.stats[[paste0("stats.", remnants[i1])]] 
            for (i2 in seq_len(i1-1L)) {
                stats2 <- ref.stats[[paste0("stats.", remnants[i2])]] 

                # Obtaining the IUT and setting opposing log-fold changes to 1.
                max.p <- pmax(stats1$p.value, stats2$p.value)
                max.p[sign(stats1$logFC) != sign(stats2$logFC)] <- 1
                    
                # Correcting across genes.
                adj.p <- p.adjust(max.p, method="BH")
                N <- sum(adj.p <= threshold, na.rm=TRUE)

                if (N < best.N) {
                    best.N <- N
                    chosen <- which.min(max.p)
                    best.gene <- chosen
                    best.p <- adj.p[chosen]
                    best.parents <- c(i1, i2)
                }
            }
        }

        print(best.N)
        print(best.parents)
        print(best.gene)
        print(best.p)
        collected[[ref]] <- DataFrame(N=best.N, 
            Source1=remnants[best.parents[1]],
            Source2=remnants[best.parents[2]],
            best=rownames(ref.stats)[best.gene], 
            p.value=best.p, row.names=ref)
    }

    # Returning the DataFrame of compiled results.
    out <- do.call(rbind, collected)
    out[order(out$N),]
}
