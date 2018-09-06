#' @importFrom scater librarySizeFactors normalize
#' @importFrom BiocGenerics "sizeFactors<-" sizeFactors
#' @importFrom stats p.adjust
#' @importFrom methods as
#' @importClassesFrom S4Vectors SimpleList
.doublet_cluster <- function(x, clusters, subset.row=NULL, threshold=0.05, ...) 
# Finds evidence that a cluster is _not_ a doublet of two other clusters.
# Absence of such evidence should be a warning flag for that cluster.
# 
# written by Aaron Lun
# created 16 April 2018
{
    if (length(unique(clusters)) < 3L) {
        stop("need at least three clusters to detect doublet clusters")
    }

    # Computing normalized counts using the library size (looking for compositional differences!)
    sce <- SingleCellExperiment(list(counts=x))
    sizeFactors(sce) <- librarySizeFactors(x, subset_row=subset.row)
    sce <- normalize(sce, return_log=TRUE)

    degs <- findMarkers(sce, clusters=clusters, subset.row=subset.row, full.stats=TRUE, ...)
    med.lib.size <- vapply(split(sizeFactors(sce), clusters), FUN=median, FUN.VALUE=0)
	n.cluster <- table(clusters)/length(clusters)

    # Setting up the output.
    all.clusters <- names(degs)
    collected.top <- collected.all <- vector("list", length(all.clusters))
    names(collected.top) <- names(collected.all) <- all.clusters

    # Running through all pairs of clusters and testing against the third cluster.
    for (ref in all.clusters) {
        ref.stats <- degs[[ref]]
        remnants <- setdiff(all.clusters, ref)

        num <- length(remnants) * (length(remnants) - 1L)/2L
        all.N <- all.gene <- all.parent1 <- all.parent2 <- integer(num)
        all.p <- numeric(num)
        idx <- 1L

        for (i1 in seq_along(remnants)) {
            stats1 <- ref.stats[[paste0("stats.", remnants[i1])]] 
            for (i2 in seq_len(i1-1L)) {
                stats2 <- ref.stats[[paste0("stats.", remnants[i2])]] 

                # Obtaining the IUT and setting opposing log-fold changes to 1.
                max.p <- pmax(stats1$p.value, stats2$p.value)
                max.p[sign(stats1$logFC) != sign(stats2$logFC)] <- 1
                    
                # Correcting across genes.
                adj.p <- p.adjust(max.p, method="BH")
                best.gene <- which.min(max.p)[1] # [1] gives NA when there are no genes, which avoids nrow mismatch in DataFrame().

                all.N[idx] <- sum(adj.p <= threshold, na.rm=TRUE)
                all.gene[idx] <- best.gene
                all.p[idx] <- adj.p[best.gene]
                all.parent1[idx] <- i1
                all.parent2[idx] <- i2
                idx <- idx + 1L
            }
        }

        # Formatting the output.
        parent1 <- remnants[all.parent1]
        parent2 <- remnants[all.parent2]

        stats <- DataFrame(source1=parent1, source2=parent2, 
            N=all.N, best=rownames(ref.stats)[all.gene], p.value=all.p,
            lib.size1=unname(med.lib.size[parent1]/med.lib.size[ref]), 
			lib.size2=unname(med.lib.size[parent2]/med.lib.size[ref]))

        o <- order(all.N)
        top <- cbind(stats[o[1],], prop=n.cluster[[ref]])
        rownames(top) <- ref
        collected.top[[ref]] <- top
        collected.all[[ref]] <- stats[o,]
    }

    # Returning the DataFrame of compiled results.
    out <- do.call(rbind, collected.top)
    out$all.pairs <- as(collected.all, "SimpleList")
    out[order(out$N),]
}

##############################
# S4 method definitions here #
##############################

#' @export
setGeneric("doubletCluster", function(x, ...) standardGeneric("doubletCluster"))

#' @export
setMethod("doubletCluster", "ANY", .doublet_cluster)

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("doubletCluster", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="counts", get.spikes=FALSE) {
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .doublet_cluster(assay(x, i=assay.type), ..., subset.row=subset.row)
})
