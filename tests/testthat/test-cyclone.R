# This tests the performance of the cyclone method, by comparing it to
# a slow R-based implementation.

# require(scran); require(testthat)

classif.single <- function(cell, markers,Nmin.couples) { 
    test<-unlist(cell[markers[,1]]-cell[markers[,2]])
    t1<-length(test[test>0])
    tot<-length(test[test!=0])
    if (tot < Nmin.couples){ return(NA) }  
    return(t1/tot)
}

random.success <- function(cell, markers, N, Nmin, Nmin.couples) {  
    cell.random <- .Call(scran:::cxx_auto_shuffle, cell, N)
    success <- apply(cell.random, 2, classif.single, markers=markers, Nmin.couples=Nmin.couples)

    success<-success[!is.na(success)]
    if(length(success)<Nmin){ return(NA) }
    
    test<-classif.single(cell,markers,Nmin.couples) 
    if(is.na(test)){ return(NA) } 
    
    return(length(success[success<test])/length(success))
}

refFUN <- function(x, pairs) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    gene.names <- rownames(x)

    chosen.x <- list()
    for (p in names(pairs)) {
        curp <- pairs[[p]]
        m1 <- match(curp$first, gene.names)
        m2 <- match(curp$second, gene.names)
        keep <- !is.na(m1) & !is.na(m2)
        m1 <- m1[keep]
        m2 <- m2[keep]
        
        all.present <- as.vector(rbind(m1, m2))
        all.present <- all.present[!duplicated(all.present)]
        chosen.x[[p]] <- x[all.present,,drop=FALSE]
        pairs[[p]] <- data.frame(first=match(m1, all.present),
                                 second=match(m2, all.present))
    }

    N <- 1000L
    Nmin <- 100L
    Nmin.couples <- 50L
    score.G1<-apply(chosen.x$G1, 2, function(x) random.success(cell=x,markers=pairs$G1,N=N,Nmin=Nmin,Nmin.couples=Nmin.couples))#, genes.list=genes.list))
    score.S<-apply(chosen.x$S, 2, function(x) random.success(cell=x,markers=pairs$S,N=N, Nmin=Nmin, Nmin.couples=Nmin.couples))#,genes.list=genes.list))
    score.G2M<-apply(chosen.x$G2M, 2, function(x) random.success(cell=x,markers=pairs$G2M,N=N, Nmin=Nmin,Nmin.couples=Nmin.couples))#,genes.list=genes.list))

    scores <- data.frame(G1=score.G1, S=score.S, G2M=score.G2M) 
    scores.normalised<-data.frame(t(apply(scores, 1, function(x) (x)/sum(x))))
    return(list(scores=scores, normalized.scores=scores.normalised))
}

####################################################################################################

# Setting cores constant to avoid messing with the seeds.

current <- bpparam()
bpworkers(current) <- 1
register(current)

# Spawning training data.

all.names <- paste0("X", seq_len(500))
Ngenes <- length(all.names)
all.pairs <- combn(Ngenes, 2)
re.pairs <- data.frame(first=all.names[all.pairs[1,]], second=all.names[all.pairs[2,]])

set.seed(100)
markers <- list(G1=re.pairs[sample(nrow(re.pairs), 100),],
                 S=re.pairs[sample(nrow(re.pairs), 200),],
               G2M=re.pairs[sample(nrow(re.pairs), 500),])

# Spawning live runs.

set.seed(1000)
Ncells <- 10
X <- matrix(rnorm(Ngenes*Ncells), ncol=Ncells)
rownames(X) <- all.names

set.seed(100)
reference <- refFUN(X, markers)

set.seed(100)
observed <- cyclone(X, markers)

expect_equal(reference$scores, observed$scores)
expect_equal(reference$normalized.scores, observed$normalized.scores)

set.seed(1001)
X <- matrix(rpois(Ngenes*Ncells, lambda=10), ncol=Ncells)
rownames(X) <- all.names

set.seed(100)
reference <- refFUN(X, markers)

set.seed(100)
observed <- cyclone(X, markers)

expect_equal(reference$scores, observed$scores)
expect_equal(reference$normalized.scores, observed$normalized.scores)

# Low counts to induce more ties.

set.seed(1002)
X <- matrix(rpois(Ngenes*Ncells, lambda=5), ncol=Ncells)
rownames(X) <- all.names

set.seed(100)
reference <- refFUN(X, markers)

set.seed(100)
observed <- cyclone(X, markers)

expect_equal(reference$scores, observed$scores)
expect_equal(reference$normalized.scores, observed$normalized.scores)

set.seed(1003)
X <- matrix(rpois(Ngenes*Ncells, lambda=1), ncol=Ncells)
rownames(X) <- all.names

set.seed(100)
reference <- refFUN(X, markers)

set.seed(100)
observed <- cyclone(X, markers)

expect_equal(reference$scores, observed$scores)
expect_equal(reference$normalized.scores, observed$normalized.scores)

# Checking that it also works with SCESet objects.

set.seed(1004)
X <- matrix(rpois(Ngenes*Ncells, lambda=100), ncol=Ncells)
rownames(X) <- all.names
X2 <- newSCESet(countData=as.data.frame(X))

set.seed(100)
reference <- refFUN(X, markers)

set.seed(100)
observed1 <- cyclone(X, markers)

set.seed(100)
observed2 <- cyclone(X2, markers, assay="counts")

expect_equal(reference$scores, observed2$scores)
expect_equal(observed1$scores, observed2$scores)
expect_equal(reference$normalized.scores, observed2$normalized.scores)
expect_equal(observed1$normalized.scores, observed2$normalized.scores)

# Odd behaviour with no cells.

out <- cyclone(X[,0], markers)
expect_identical(nrow(out$scores), 0L)
expect_identical(colnames(out$scores), c("G1", "S", "G2M"))
expect_identical(nrow(out$normalized.scores), 0L)
expect_identical(colnames(out$normalized.scores), c("G1", "S", "G2M"))

# Odd behaviour with no markers.

no.markers <- list(G1=re.pairs[0,],
                    S=re.pairs[0,],
                  G2M=re.pairs[0,])
out <- cyclone(X, no.markers)
expect_identical(colnames(out$scores), c("G1", "S", "G2M"))
expect_identical(nrow(out$scores), ncol(X))
expect_true(all(is.na(out$scores)))
expect_identical(colnames(out$normalized.scores), c("G1", "S", "G2M"))
expect_identical(nrow(out$normalized.scores), ncol(X))
expect_true(all(is.na(out$normalized.scores)))


