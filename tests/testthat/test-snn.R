# Checks the construction of the SNN graph.
# require(scran); require(testthat); source("test-snn.R")

# Constructing a reference value.

library(igraph)
library(FNN)
check <- function(vals, k=10) {
    g <- buildSNNGraph(vals, k=k)
    nn.out <- get.knn(t(vals), k=k)
    IDX <- cbind(seq_len(ncol(vals)), nn.out$nn.index)

    ncells <- ncol(vals)
    expect_identical(seq_len(ncells), as.vector(V(g)))
    for (i in seq_len(ncells)) { 
        inn <- IDX[i,]
        collected <- numeric(ncells)

        for (j in seq_len(ncells)) {
            jnn <- IDX[j,]
            shared <- intersect(inn, jnn)
            if (length(shared)==0) next
            s <- k + 1 - 0.5*(match(shared, inn) + match(shared, jnn))
            collected[j] <- max(s)
        }
        collected[i] <- 0
        expect_equal(collected, g[i])
    }
    return(NULL)
}

set.seed(20000)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)

dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
check(dummy, k=10)

dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
check(dummy, k=20)

dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
check(dummy, k=5)

# Checking that the value is sensible with subset.row.

are_graphs_same <- function(g1, g2) {
    expect_equal(g1[], g2[])
    return(TRUE)
}

selected <- sample(ngenes, 50)
g <- buildSNNGraph(dummy[selected,])
g2 <- buildSNNGraph(dummy, subset.row=selected)
are_graphs_same(g, g2)

# Checking SCESet construction.

suppressWarnings(sce <- newSCESet(countData=2^dummy))
g <- buildSNNGraph(sce)
g2 <- buildSNNGraph(exprs(sce))
are_graphs_same(g, g2)

g <- buildSNNGraph(sce, assay="counts")
g2 <- buildSNNGraph(2^dummy)
are_graphs_same(g, g2)

g <- buildSNNGraph(sce, subset.row=selected)
g2 <- buildSNNGraph(sce[selected,])
are_graphs_same(g, g2)

sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=selected))
setSpike(sce) <- "ERCC"
g <- buildSNNGraph(sce)
g2 <- buildSNNGraph(sce[-selected,])
are_graphs_same(g, g2)
