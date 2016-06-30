find.markers <- function(data1, data2, data3, gene.names, fraction=0.5)
# This identifies pairs of genes whose relative expression is > 0 in 
# at least a 'fraction' of cells in one phase is < 0 in at least 
# 'fraction' of the cells in each of the other phases.
{
    Ngenes <- ncol(data1)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' vector must be equal to 'x' nrows")
    }
    if (Ngenes!=ncol(data2) || Ngenes!=ncol(data3)) { 
        stop("number of genes in each phase must be the same")
    }
    if (nrow(data1)==0L || nrow(data2)==0L || nrow(data3)==0L) {
        stop("each phase must have at least one cell")
    }

    # Calculating thresholds.
    Nthr1 <- ceiling(nrow(data1) * fraction)
    Nthr2 <- ceiling(nrow(data2) * fraction)
    Nthr3 <- ceiling(nrow(data3) * fraction)

    if (Ngenes) { 
        collected <- list()
        counter <- 1L
        for (i in seq_len(Ngenes-1L)) { 
            others <- (i+1):Ngenes
            diff1 <- data1[,i] - data1[,others,drop=FALSE]
            diff2 <- data2[,i] - data2[,others,drop=FALSE]
            diff3 <- data3[,i] - data3[,others,drop=FALSE]

            npos1 <- colSums(diff1 > 0)
            npos2 <- colSums(diff2 > 0)
            npos3 <- colSums(diff3 > 0)
            nneg1 <- colSums(diff1 < 0)
            nneg2 <- colSums(diff2 < 0)
            nneg3 <- colSums(diff3 < 0)

            chosen <- others[npos1 >= Nthr1 & nneg2 >= Nthr2 & nneg3 >= Nthr3]
            if (length(chosen)) { 
                collected[[counter]] <- cbind(i, chosen)
                counter <- counter + 1L
            }
            chosen.flip <- others[nneg1 >= Nthr1 & npos2 >= Nthr2 & npos3 >= Nthr3]
            if (length(chosen.flip)) { 
                collected[[counter]] <- cbind(chosen.flip, i)
                counter <- counter + 1L
            }
        }
        collected <- do.call(rbind, collected)
        g1 <- gene.names[collected[,1]]
        g2 <- gene.names[collected[,2]]
    } else {
        g1 <- g2 <- character(0)
    }

    return(data.frame(first=g1, second=g2, stringsAsFactors=FALSE))
}

setGeneric("sandbag", function(x, ...) standardGeneric("sandbag"))

setMethod("sandbag", "matrix", function(x, is.G1, is.S, is.G2M, gene.names=rownames(x), fraction=0.5, subset.row=NULL) 
# Identifies the relevant pairs before running 'cyclone'.
# Basically runs through all combinations of 'find.markers' for each phase. 
#
# written by Aaron Lun
# based on code by Antonio Scialdone
# created 22 January 2016 
# last modified 8 June 2016
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    data.G1 <- t(x[subset.row,is.G1,drop=FALSE])
    data.S <- t(x[subset.row,is.S,drop=FALSE])
    data.G2M <- t(x[subset.row,is.G2M,drop=FALSE]) 
    gene.names <- gene.names[subset.row]

    G1.marker.pairs <- find.markers(data.G1, data.S, data.G2M, fraction=fraction, gene.names=gene.names)
    S.marker.pairs <- find.markers(data.S, data.G1, data.G2M, fraction=fraction, gene.names=gene.names)
    G2M.marker.pairs <- find.markers(data.G2M, data.G1, data.S, fraction=fraction, gene.names=gene.names)
    return(list(G1=G1.marker.pairs, S=S.marker.pairs, G2M=G2M.marker.pairs))
})

setMethod("sandbag", "SCESet", function(x, is.G1, is.S, is.G2M, subset.row=NULL, ..., assay="counts", get.spikes=FALSE) {
    if (is.null(subset.row)) { 
        subset.row <- .spikeSubset(x, get.spikes)
    }
    sandbag(assayDataElement(x, assay), is.G1=is.G1, is.S=is.S, is.G2M=is.G2M, ..., subset.row=subset.row)
})
