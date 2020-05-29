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
    centers <- rbind(A=c(1, 0), B=c(3, 0), C=c(4.5, 0), D=c(6, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y))

    mst <- createClusterMST(centers)
    ordering <- orderClusterMST(y, clusters, centers, mst)

    m <- match(clusters, rownames(centers))
    left.lim <- c(min(centers[,1]), centers[,1])[m]
    right.lim <- c(centers[,1], max(centers[,1]))[m+1]
    expect_equivalent(as.numeric(ordering), pmin(right.lim, pmax(left.lim, y[,1])) - 1)
})

set.seed(99000002)
test_that("MST ordering works as expected for branching", {
    centers <- rbind(base=c(0, 0), A=c(1, 0), B=c(0, 2), C=c(0, -3), D=c(-4, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + runif(length(y), -0.4, 0.4)

    mst <- createClusterMST(centers)
    ordering <- orderClusterMST(y, clusters, centers, mst, start="base")

    # Everyone gets assigned to only one pseudotime.
    expect_identical(ncol(ordering), nrow(centers) - 1L)
    expect_true(all(rowSums(!is.na(ordering))==1))

    # Computing the branched pseudo-times.
    expect_equivalent(ordering[clusters=="A",1], pmin(y[clusters=="A",1], 1))
    expect_equivalent(ordering[clusters=="B",2], pmin(y[clusters=="B",2], 2))
    expect_equivalent(ordering[clusters=="C",3], -pmax(y[clusters=="C",2], -3))
    expect_equivalent(ordering[clusters=="D",4], -pmax(y[clusters=="D",1], -4))
})

set.seed(99000002)
test_that("MST ordering works as expected for a complex case", {
    centers <- rbind(base=c(0, 0), A=c(1, 0), B=c(0, 2), C=c(0, -3), D=c(-4, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + runif(length(y), -0.4, 0.4)

    mst <- createClusterMST(centers)
    ordering <- orderClusterMST(y, clusters, centers, mst, start="base")

    # Everyone gets assigned to only one pseudotime.
    expect_identical(ncol(ordering), nrow(centers) - 1L)
    expect_true(all(rowSums(!is.na(ordering))==1))

    # Computing the branched pseudo-times.
    expect_equivalent(ordering[clusters=="A",1], pmin(y[clusters=="A",1], 1))
    expect_equivalent(ordering[clusters=="B",2], pmin(y[clusters=="B",2], 2))
    expect_equivalent(ordering[clusters=="C",3], -pmax(y[clusters=="C",2], -3))
    expect_equivalent(ordering[clusters=="D",4], -pmax(y[clusters=="D",1], -4))

    # Computing the expected partition for the base.
    base <- clusters=="base"
    expect_equivalent(
        rowMeans(ordering[base,,drop=FALSE], na.rm=TRUE),
        pmax(abs(y[base,1]), abs(y[base,2]))
    )
})

set.seed(990000021)
test_that("MST ordering shares values at zero", {
    X <- as.matrix(expand.grid(-1:1, -1:1))
    Y <- rbind(A=c(-1, -1), B=c(0, 0), C=c(1, -1))
    mst <- createClusterMST(Y)
    ordering <- orderClusterMST(X, rep("B", nrow(X)), Y, mst, start="B")

    expect_identical(sum(ordering==0, na.rm=TRUE)/ncol(ordering), 4)
    expect_identical(which(ordering[,1]==0), which(ordering[,2]==0))

    # Still the case when it's not happening on the starting node.
    Y2 <- rbind(Y, D=c(0, -1))
    mst <- igraph::make_graph(c("D", "B", "B", "A", "B", "C"))
    ordering2 <- orderClusterMST(X, rep("B", nrow(X)), Y2, mst, start="D")

    expect_identical(sum(ordering2==0, na.rm=TRUE)/ncol(ordering2), 1)
    expect_identical(sum(ordering2==1, na.rm=TRUE)/ncol(ordering2), 4)

    expect_identical(which(ordering2[,1]==1), which(ordering2[,2]==1))
    expect_identical(which(ordering2[,1]==0), which(ordering2[,2]==0))
})

set.seed(99000003)
test_that("MST ordering is robust to transformations", {
    centers <- matrix(rnorm(20, sd=3), ncol=2)
    rownames(centers) <- letters[1:10]

    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y), 0.1)

    mst <- createClusterMST(centers)
    ref <- orderClusterMST(y, clusters, centers, mst)

    # For your viewing pleasure:
    # plot(y[,1], y[,2], col=topo.colors(21)[cut(rowMeans(ref, na.rm=TRUE), 21)])
    # stuff <- connectClusterMST(centers, mst, combined=FALSE)
    # segments(stuff$start$dim1, stuff$start$dim2, stuff$end$dim1, stuff$end$dim2, lwd=5)

    # Applying various transformation.
    expect_equivalent(ref*2, orderClusterMST(y*2, clusters, centers*2, mst))
    expect_equivalent(ref, orderClusterMST(y+1, clusters, centers+1, mst))

    rotation <- matrix(c(cos(pi/4), sin(pi/4), -sin(pi/4), cos(pi/4)), ncol=2)
    expect_equivalent(ref, orderClusterMST(y %*% rotation, clusters, centers %*% rotation, mst))
})
