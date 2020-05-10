# Checks the construction of the SNN graph.
# require(scran); require(testthat); source("setup.R"); source("test-build-snn.R")

library(igraph)
check <- function(vals, k=10, type="rank")
# Checking against a slow reference calculator.
{
    g <- buildSNNGraph(vals, k=k, d=NA, type=type) # turning off PCA.
    nn.out <- BiocNeighbors::findKNN(t(vals), k=k)
    IDX <- cbind(seq_len(ncol(vals)), nn.out$index)

    ncells <- ncol(vals)
    expect_identical(seq_len(ncells), as.vector(V(g)))
    for (i in seq_len(ncells)) { 
        inn <- IDX[i,]
        collected <- numeric(ncells)

        for (j in seq_len(ncells)) {
            jnn <- IDX[j,]
            shared <- intersect(inn, jnn)
            if (length(shared)==0) {
                next
            }
            if (type=="rank") {
                s <- k + 1 - 0.5*(match(shared, inn) + match(shared, jnn))
                collected[j] <- max(s)
            } else if (type=="number") {
                collected[j] <- length(shared)
            } else {
                collected[j] <- length(shared) / length(union(inn, jnn))
            }
        }
        collected[i] <- 0
        expect_equal(collected, g[i])
    }
    return(NULL)
}

set.seed(20000)
ncells <- 200
ngenes <- 500

test_that("buildSNNGraph gives same results as a reference", {
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=10)
    
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=20)
    
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=5)

    # Checking 'number' mode.  
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=10, type="number")
    
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=20, type="number")
    
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=5, type="number")

    # Checking 'jaccard' mode.  
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=10, type="jaccard")
    
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=20, type="jaccard")
    
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    check(dummy, k=5, type="jaccard")
})

# Checking that the value is sensible with subset.row.

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

# Checking SCESet construction.

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

# Checking PCA was working.

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

    # Unaffected by spike-in and subset.row specifications (correctly).
    sub <- sample(ngenes, 50)
    alt <- buildSNNGraph(X, use.dimred="PCA", subset.row=sub)
    are_graphs_same(ref, alt)
})

# Silly inputs.

test_that("buildSNNGraph fails on silly inputs", {
    dummy <- matrix(rnorm(ngenes*20), ncol=20, nrow=ngenes)
    expect_warning(out <- buildSNNGraph(dummy, k=50, d=NA), "capped")
    expect_warning(out2 <- buildSNNGraph(dummy, k=ncol(dummy)-1L, d=NA), NA)
    are_graphs_same(out, out2)

    expect_error(buildSNNGraph(dummy[0,], d=NA), NA) # shouldn't fail, but shouldn't generate anything particularly useful.

    expect_warning(expect_error(buildSNNGraph(dummy[,0], d=NA), "must be positive"), "capped")
})

# Checking that buildKNNGraph also works.

KMAKE <- function(dummy, k, directed=FALSE) { 
    ncells <- ncol(dummy)
    collated <- matrix(0, ncells, ncells)
    for (cell in seq_len(ncells)) {
        d2 <- colSums((dummy[,cell] - dummy)^2)
        chosen <- setdiff(order(d2), cell)[seq_len(k)]
        collated[cell,chosen] <- 1
        if (!directed) { 
            collated[chosen,cell] <- 1
        }
    }
    return(as(collated, "dgCMatrix"))
}

test_that("buildKNNGraph works correctly", {
    dummy <- matrix(rnorm(ngenes*ncells), ncol=ncells, nrow=ngenes)
    g <- buildKNNGraph(dummy, k=10, d=NA)
    expect_false(is.directed(g))
    expect_equal(g[], KMAKE(dummy, k=10))

    g <- buildKNNGraph(dummy, k=10, d=NA, directed=TRUE)
    expect_true(is.directed(g))
    expect_equal(g[], KMAKE(dummy, k=10, directed=TRUE))

    g <- buildKNNGraph(dummy, k=20, d=NA)
    expect_equal(g[], KMAKE(dummy, k=20))

    pc.out <- prcomp(t(dummy), rank.=10)
    g <- buildKNNGraph(dummy, k=20, d=10)
    expect_equal(g[], KMAKE(t(pc.out$x), k=20))
})

# Avoid normalize() overwriting scuttle's normalize() in other files.

detach("package:igraph", character.only=TRUE)

