find.markers <- function(id1, id2, id3, training.data, fraction=0.5, gene.names=rownames(training.data))
# This identifies pairs of genes whose relative expression is > 0 in 
# at least a 'fraction' of cells in one phase is < 0 in at least 
# 'fraction' of the cells in each of the other phases.
{
    Ngenes <- nrow(training.data)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' vector must be equal to 'training.data' nrows")
    }

    data1 <- t(training.data[,id1,drop=FALSE])
    data2 <- t(training.data[,id2,drop=FALSE])
    data3 <- t(training.data[,id3,drop=FALSE]) 
    if (nrow(data1)==0L || nrow(data2)==0L || nrow(data3)==0L) {
        stop("each phase must have at least one cell")
    }

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

classify.single <- function(cell, markers, Nmin.couples) 
# Count number of successes for a single set of markers, given a cell's expression values.
# Successes are defined as those pairs where the first gene is more highly expressed.
{
    test <- cell[markers[,1]] - cell[markers[,2]]
    t1 <- sum(test>0)
    tot <- sum(test!=0)
    if (tot < Nmin.couples) { return(NA) }
    return(t1/tot)
}

random.success <- function(cell, markers, N, Nmin, Nmin.couples)
# Given a set of markers, a number N of trials and a cell, find the number of hits in each of the N randomised set of markers.
# This returns the probability of randoml obtaining a fraction of hits lower than the observed fraction.
# The null hypothesis is that the expression of each gene is sampled from the empirical distribution in 'cell',
# in a manner that is independent of the pairings between genes. We then calculate the classification based on those pairings.
{
    success <- sapply(seq_len(N), function(x) {    
                      cell.random <- cell[sample(length(cell))]
                      classify.single(cell.random, markers, Nmin.couples)
    })
  
    success <- success[!is.na(success)]
    test <- classify.single(cell, markers, Nmin.couples)

    if(length(success) < Nmin || is.na(test)) { 
        warning("not enough gene pairs with different expression values")
        return(NA) 
    }
    return(sum(success<test)/length(success))
}

setGeneric("sandbag", function(x, ...) standardGeneric("sandbag"))

setMethod("sandbag", "ANY", function(x, is.G1, is.S, is.G2M, gene.names=rownames(x), fraction=0.5) 
# Identifies the relevant pairs before running 'cyclone'.
# Basically runs through all combinations of 'find.markers' for each phase. 
#
# written by Aaron Lun
# based on code by Antonio Scialdone
# created 22 January 2016 
# last modified 17 February 2016          
{
    x <- as.matrix(x)
    G1.marker.pairs <- find.markers(id1=is.G1, id2=is.S, id3=is.G2M, training.data=x, fraction=fraction, gene.names=gene.names)
    S.marker.pairs <- find.markers(id1=is.S, id2=is.G1, id3=is.G2M, training.data=x, fraction=fraction, gene.names=gene.names)
    G2M.marker.pairs <- find.markers(id1=is.G2M, id2=is.G1, id3=is.S, training.data=x, fraction=fraction, gene.names=gene.names)
    return(list(G1=G1.marker.pairs, S=S.marker.pairs, G2M=G2M.marker.pairs))
})

setMethod("sandbag", "SCESet", function(x, ..., assay="counts", get.spikes=FALSE) {
    sandbag(.getUsedMatrix(x, assay, get.spikes), ...)
})
