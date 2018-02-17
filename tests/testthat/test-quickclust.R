# This tests out the quickCluster method in scran.
# require(scran); require(testthat); source("test-quickclust.R")

set.seed(30000)
ncells <- 700
ngenes <- 1000

dummy <- matrix(rnbinom(ncells*ngenes, mu=10, size=20), ncol=ncells, nrow=ngenes)
known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0
    
test_that("quickCluster works in the simple case", {
    emp.clusters <- quickCluster(dummy)
    expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
    shuffled <- c(1:50, 301:350, 601:650)
    expect_identical(quickCluster(dummy, subset.row=shuffled), emp.clusters)
})

test_that("quickCluster correctly computes the ranks", {
    emp.ranks <- quickCluster(dummy, get.ranks=TRUE)
    ref <- apply(dummy, 2, FUN=function(y) {
        r <- rank(y)
        r <- r - mean(r)
        r/sqrt(sum(r^2))/2
    })
    expect_equal(emp.ranks, ref)
    
    shuffled <- c(50:100, 401:350, 750:850)
    emp.ranks <- quickCluster(dummy, get.ranks=TRUE, subset.row=shuffled)
    ref <- apply(dummy, 2, FUN=function(y) {
        r <- rank(y[shuffled])
        r <- r - mean(r)
        r/sqrt(sum(r^2))/2
    })
    expect_equal(emp.ranks, ref)
})

# Checking out that clustering is consistent with that based on correlations.

set.seed(300001)
test_that("quickCluster is consistent with clustering on correlations", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    obs <- quickCluster(mat)
    
    refM <- sqrt(0.5*(1 - cor(mat, method="spearman")))
    distM <- as.dist(refM) 
    htree <- hclust(distM, method='ward.D2')
    clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=200, distM=refM, verbose=0))
    expect_identical(clusters, as.integer(obs)) # this can complain if unassigned, as 0 becomes 1 in as.integer().

    # Works with min.size set.
    obs <- quickCluster(mat, min.size=50)
    clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=refM, verbose=0))
    expect_identical(clusters, as.integer(obs))

    # Works with max.size set.
    obs <- quickCluster(mat, min.size=20)
    obs2 <- quickCluster(mat, min.size=20, max.size=40)
    expect_false(all(table(obs) <= 40))
    expect_true(all(table(obs2) <= 40))
    expect_identical(length(unique(paste(obs, obs2, sep="."))), nlevels(obs2)) # obs2 is fully nested in obs.
})

test_that("quickCluster functions correctly with subsetting", {
    # Works properly with subsetting.
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    subset.row <- 15:1 
    obs <- quickCluster(mat, min.size=50, subset.row=subset.row)
    expect_identical(obs, quickCluster(mat[subset.row,], min.size=50)) # Checking that it behaves properly.
    expect_false(identical(quickCluster(mat), obs)) # It should return different results.

    refM <- sqrt(0.5*(1 - cor(mat[subset.row,], method="spearman")))
    distM <- as.dist(refM) 
    htree <- hclust(distM, method='ward.D2')
    clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=refM, verbose=0))
    expect_identical(clusters, as.integer(obs))

    # Handles the mean.
    obs <- quickCluster(mat, min.size=50, min.mean=5)
    expect_identical(obs, quickCluster(mat, min.size=50, subset.row=scater::calcAverage(mat) >= 5))
    rnks <- quickCluster(mat, get.ranks=TRUE, min.mean=5)
    expect_identical(rnks, quickCluster(mat, get.ranks=TRUE, subset.row=scater::calcAverage(mat) >= 5))
})

# Checking that we're executing the igraph methods correctly.

set.seed(300002)
test_that("quickCluster with igraph works with min.size settings", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    obs <- quickCluster(mat, min.size=0, method="igraph")
    expect_identical(length(obs), ncol(mat))

    ref <- quickCluster(mat, get.rank=TRUE)
    snn <- buildSNNGraph(ref, pc.approx=TRUE) 
    library(igraph)
    out <- cluster_fast_greedy(snn)
    expect_identical(factor(out$membership), obs)
    
    min.size <- 100 # Checking that min.size merging works.
    expect_false(all(table(obs) >= min.size))
    obs2 <- quickCluster(mat, min.size=min.size, method="igraph")
    expect_true(all(table(obs2) >= min.size))

    combined <- paste0(obs, ".", obs2)
    expect_identical(length(unique(combined)), length(unique(obs))) # Confirm that they are nested.

    # Checking that arguments are passed along.
    obs <- quickCluster(mat, min.size=0, method="igraph", pc.approx=FALSE)
    snn <- buildSNNGraph(ref, pc.approx=FALSE) 
    out <- cluster_fast_greedy(snn)
    expect_identical(factor(out$membership), obs)

    obs <- quickCluster(mat, min.size=0, method="igraph", d=NA)
    snn <- buildSNNGraph(ref, d=NA)
    out <- cluster_fast_greedy(snn)
    expect_identical(factor(out$membership), obs)
})


# Creating an example where quickCluster should behave correctly.

ncells <- 600
ngenes <- 200
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- seq_len(ngenes)/100
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:40,known.clusters==1L] <- 0
dummy[41:80,known.clusters==2L] <- 0  
dummy[81:120,known.clusters==3L] <- 0

test_that("quickCluster reports the correct clusters", {
    out <- quickCluster(dummy)          
    expect_identical(length(unique(paste(out, known.clusters))), 3L)
}) 

test_that("quickCluster warns or fails on silly inputs", {
    expect_error(quickCluster(dummy[0,]), "rank variances of zero detected for a cell")
    expect_error(quickCluster(dummy[,0]), "fewer cells than the minimum cluster size")

    # Checking for a warning upon unassigned cells.    
    leftovers <- 100
    expect_warning(forced <- quickCluster(dummy[,c(which(known.clusters==1), 
                                                   which(known.clusters==2), 
                                                   which(known.clusters==3)[seq_len(leftovers)])]), 
                   sprintf("%i cells were not assigned to any cluster", leftovers))
    expect_identical(as.character(tail(forced, leftovers)), rep("0", leftovers))
})

set.seed(20002)
test_that("quickCluster works on SingleCellExperiment objects", {
    dummy <- matrix(rpois(50000, lambda=5), nrow=50)
    rownames(dummy) <- paste0("X", seq_len(nrow(dummy)))
    X <- SingleCellExperiment(list(counts=dummy))
    emp.clusters <- quickCluster(X)
    expect_identical(emp.clusters, quickCluster(counts(X)))

    # Checking correct interplay between spike-ins and subset.row.
    isSpike(X, "ERCC") <- 1:20
    expect_identical(quickCluster(X), quickCluster(counts(X)[-(1:20),]))
    subset.row <- 1:25*2
    expect_identical(quickCluster(X, subset.row=subset.row), 
                     quickCluster(counts(X)[setdiff(subset.row, 1:20),]))
    expect_identical(quickCluster(X, subset.row=subset.row, get.spikes=TRUE), 
                     quickCluster(counts(X)[subset.row,]))   
})
