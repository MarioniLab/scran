# Checks the construction of the SNN graph.
# require(scran); require(testthat); source("setup.R"); source("test-build-snn.R")

ngenes <- 500
ncells <- 200

are_graphs_same <- function(g1, g2) {
    expect_equal(g1[], g2[])
    return(TRUE)
}

set.seed(20001)
test_that("Subsetting does not change the result", {
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)

    selected <- sample(ngenes, 50)
    g <- buildSNNGraph(dummy[selected,])
    g2 <- buildSNNGraph(dummy, subset.row=selected)
    are_graphs_same(g, g2)

    selected <- rbinom(ngenes, 1, 0.5)==1
    g <- buildSNNGraph(dummy[selected,])
    g2 <- buildSNNGraph(dummy, subset.row=selected)
    are_graphs_same(g, g2)
})

set.seed(20002)
test_that("buildSNNGraph works properly on SingleCellExperiment objects", {
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    sce <- SingleCellExperiment(list(counts=2^dummy, logcounts=dummy))
    g <- buildSNNGraph(sce)
    g2 <- buildSNNGraph(assay(sce, "logcounts"))
    are_graphs_same(g, g2)
    
    g <- buildSNNGraph(sce, assay.type="counts")
    g2 <- buildSNNGraph(assay(sce, "counts"))
    are_graphs_same(g, g2)
    
    selected <- sample(ngenes, 50)
    g <- buildSNNGraph(sce, subset.row=selected)
    g2 <- buildSNNGraph(sce[selected,])
    are_graphs_same(g, g2)
})

set.seed(20004)
test_that("buildSNNGraph with PCA works correctly", {
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    pc <- prcomp(t(dummy))
    ref <- buildSNNGraph(t(pc$x[,1:20]), k=10, d=NA)
    alt <- buildSNNGraph(dummy, k=10, d=20)
    are_graphs_same(ref, alt)
    
    ref <- buildSNNGraph(t(pc$x[,1:50]), k=10, d=NA)
    alt <- buildSNNGraph(dummy, k=10, d=50)
    are_graphs_same(ref, alt)

    # Checking that it correctly extracts stuff from the reducedDimension slot.
    X <- SingleCellExperiment(list(logcounts=dummy))
    reducedDim(X, "PCA") <- pc$x[,1:50]
    alt <- buildSNNGraph(X, use.dimred="PCA")
    are_graphs_same(ref, alt)

    # Unaffected by subset.row specifications (correctly).
    sub <- sample(ngenes, 50)
    alt <- buildSNNGraph(X, use.dimred="PCA", subset.row=sub)
    are_graphs_same(ref, alt)
})
