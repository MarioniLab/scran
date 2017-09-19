# Checks the construction of the SNN graph.
# require(scran); require(testthat); source("test-snn.R")

# Constructing a reference value.

library(igraph)
library(FNN)
check <- function(vals, k=10) {
    g <- buildSNNGraph(vals, k=k, d=NA) # turning off PCA.
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

test_that("buildSNNGraph works properly on SingleCellExperiment objects", {
    sce <- SingleCellExperiment(list(counts=2^dummy, logcounts=dummy))
    g <- buildSNNGraph(sce)
    g2 <- buildSNNGraph(assay(sce, "logcounts"))
    are_graphs_same(g, g2)
    
    g <- buildSNNGraph(sce, assay.type="counts")
    g2 <- buildSNNGraph(assay(sce, "counts"))
    are_graphs_same(g, g2)
    
    g <- buildSNNGraph(sce, subset.row=selected)
    g2 <- buildSNNGraph(sce[selected,])
    are_graphs_same(g, g2)
    
    isSpike(sce, "ERCC") <- selected
    g <- buildSNNGraph(sce)
    g2 <- buildSNNGraph(sce[-selected,])
    are_graphs_same(g, g2)
})

# Checking multi-core processing works.

dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
ref <- scran:::.find_knn(dummy, k=10, BPPARAM=SerialParam())
alt <- scran:::.find_knn(dummy, k=10, BPPARAM=SerialParam(), force=TRUE)
expect_equal(ref, alt)

dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
ref <- scran:::.find_knn(dummy, k=5, BPPARAM=SerialParam())
alt <- scran:::.find_knn(dummy, k=5, BPPARAM=SerialParam(), force=TRUE)
expect_equal(ref, alt)

# ref <- scran:::.find_knn(dummy, k=5)
# alt <- scran:::.find_knn(dummy, k=5, BPPARAM=MulticoreParam(3))
# expect_equal(ref, alt)

# Checking PCA was working.

test_that("buildSNNGRaph with PCA works correctly", {
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    pc <- prcomp(t(dummy))
    ref <- buildSNNGraph(t(pc$x[,1:20]), k=10, d=NA)
    alt <- buildSNNGraph(dummy, k=10, d=20)
    are_graphs_same(ref, alt)
    
    ref <- buildSNNGraph(t(pc$x[,1:50]), k=10, d=NA)
    alt <- buildSNNGraph(dummy, k=10, d=50)
    are_graphs_same(ref, alt)

    # Testing with IRLBA.
    set.seed(100)
    ipc <- irlba::prcomp_irlba(t(dummy - rowMeans(dummy)), n=10, center=FALSE)
    refi <- buildSNNGraph(t(ipc$x), k=10, d=NA)
    alti <- buildSNNGraph(t(ipc$x), k=10, d=10, pc.approx=TRUE, rand.seed=100)
    are_graphs_same(refi, alti)

    set.seed(200)
    ipc <- irlba::prcomp_irlba(t(dummy - rowMeans(dummy)), n=20, center=FALSE)
    refi <- buildSNNGraph(t(ipc$x), k=10, d=NA)
    alti <- buildSNNGraph(t(ipc$x), k=10, d=20, pc.approx=TRUE, rand.seed=200)
    are_graphs_same(refi, alti)

    # Checking that it correctly extracts stuff from the reducedDimension slot.
    X <- SingleCellExperiment(list(logcounts=dummy))
    reducedDim(X, "PCA") <- pc$x[,1:50]
    alt <- buildSNNGraph(X, use.dimred="PCA")
    are_graphs_same(ref, alt)

    # Ignores spike-in and subset.row specifications (correctly).
    isSpike(X, "ERCC") <- selected
    alt <- buildSNNGraph(X, use.dimred="PCA")
    are_graphs_same(ref, alt)

})

# Silly inputs.

dummy <- matrix(rnorm(ngenes*20), ncol=20, nrow=ngenes)
suppressWarnings(expect_error(buildSNNGraph(dummy[,0], d=NA), "cannot create empty graph with negative number of vertices"))
suppressWarnings(expect_error(buildSNNGraph(dummy[0,], d=NA), "cannot create empty graph with negative number of vertices"))
expect_warning(out <- buildSNNGraph(dummy, k=50, d=NA), "'k' set to the number of cells minus 1")
are_graphs_same(out, buildSNNGraph(dummy, k=19))

# Avoid normalize() overwriting scater's normalize() in other files.

detach("package:igraph", character.only=TRUE)

