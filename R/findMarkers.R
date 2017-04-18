.findMarkers <- function(x, clusters, design=NULL, subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors in 'design'.
#
# written by Aaron Lun
# created 22 March 2017
# last modified 17 April 2017    
{
    # Creating a design matrix.
    clusters <- as.factor(clusters)
    full.design <- model.matrix(~0 + clusters)
    colnames(full.design) <- clust.vals <- levels(clusters)

    if (!is.null(design)) {
        # Removing terms to avoid linearly dependencies on the intercept.
        out <- qr.solve(design, cbind(rep(1, nrow(design))))
        to.drop <- abs(out) > 1e-8
        if (any(to.drop)) {
            design <- design[,-which(to.drop)[1],drop=FALSE]
        }
        full.design <- cbind(full.design, design) # Other linear dependencies will trigger warnings.
    }
    
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    lfit <- lmFit(x[subset.row,,drop=FALSE], full.design)
    output <- vector("list", length(clust.vals))
    names(output) <- clust.vals  

    for (host in clust.vals) { 
        not.host <- clust.vals!=host
        targets <- clust.vals[not.host]
        all.p <- all.lfc <- vector("list", length(targets))
        names(all.p) <- names(all.lfc) <- targets
       
        con <- matrix(0, ncol(full.design), length(clust.vals))
        diag(con) <- -1
        con[which(!not.host),] <- 1
        con <- con[,not.host,drop=FALSE]
        colnames(con) <- targets

        fit2 <- contrasts.fit(lfit, con)
        fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
        
        for (target in targets) { 
            res <- topTable(fit2, number=Inf, sort.by="none", coef=target)
            all.p[[target]] <- res$P.Value
            all.lfc[[target]] <- res$logFC
        }

        # Computing Simes' p-value in a fully vectorised manner.
        com.p <- do.call(rbind, all.p)
        ncon <- nrow(com.p)
        ngenes <- ncol(com.p)
        gene.id <- rep(seq_len(ngenes), each=ncon)
        penalty <- rep(ncon/seq_len(ncon), ngenes) 
        o <- order(gene.id, com.p)
        com.p[] <- com.p[o]*penalty
        com.p <- t(com.p)
        smallest <- (max.col(-com.p) - 1) * ngenes + seq_len(ngenes)
        adj.min.p <- com.p[smallest]

        collected.ranks <- lapply(all.p, rank, ties="first")
        min.rank <- do.call(pmin, collected.ranks)
        marker.set <- data.frame(Top=min.rank, Gene=rownames(x)[subset.row], 
                                 FDR=p.adjust(adj.min.p, method="BH"), do.call(cbind, all.lfc), 
                                 stringsAsFactors=FALSE, check.names=FALSE)
        marker.set <- marker.set[order(marker.set$Top),]
        rownames(marker.set) <- NULL
        output[[host]] <- marker.set
    }

    return(output)
}

setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

setMethod("findMarkers", "matrix", .findMarkers)

setMethod("findMarkers", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { subset.row <- .spikeSubset(x, get.spikes) }
    .findMarkers(assayDataElement(x, assay), ..., subset.row=subset.row)
})                                 


