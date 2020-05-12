# This tests the suite of *ClusterMST functions.
# library(testthat); library(scran); source("test-cluster-mst.R")

test_that("MST construction works as expected", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    mst <- createClusterMST(y)

    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("A", "D"))
    expect_identical(vertices[igraph::degree(mst)==2], c("B", "C"))

    y <- rbind(A=c(0, 1), B=c(0, 0), C=c(1, 1), D=c(-1, 1)) 
    mst <- createClusterMST(y)

    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C", "D"))
    expect_identical(vertices[igraph::degree(mst)==3], "A")
})

test_that("MST segment reporting works as expected", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    mst <- createClusterMST(y)

    out <- connectClusterMST(y, mst)
    expect_identical(colnames(out), c("edge", "dim1", "dim2"))
    expect_identical(nrow(out), length(igraph::E(mst)) * 2L)

    named <- y
    colnames(named) <- c("PC1", "PC2")
    out2 <- connectClusterMST(named, mst)
    expect_identical(colnames(out2), c("edge", "PC1", "PC2"))

    out <- connectClusterMST(y, mst, combined=FALSE)
    expect_identical(colnames(out$start), c("edge", "dim1", "dim2"))
    expect_identical(colnames(out$end), c("edge", "dim1", "dim2"))
    expect_identical(nrow(out$start), length(igraph::E(mst)))
    expect_identical(nrow(out$end), length(igraph::E(mst)))
})

set.seed(99000001)
test_that("MST ordering works as expected for straight lines", {
    centers <- rbind(A=c(1, 0), B=c(2, 0), C=c(3, 0), D=c(4, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y))

    mst <- createClusterMST(centers)
    ordering <- orderClusterMST(y, clusters, centers, mst)

    left.lim <- pmax(centers[clusters,1] - 1, min(centers[,1]))
    right.lim <- pmin(centers[clusters,1] + 1, max(centers[,1]))
    expect_equivalent(as.numeric(ordering), pmin(right.lim, pmax(left.lim, y[,1])) - 1)
})
