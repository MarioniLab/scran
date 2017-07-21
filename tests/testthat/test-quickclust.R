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
mat <- matrix(rpois(10000, lambda=5), nrow=20)
obs <- quickCluster(mat)

refM <- sqrt(0.5*(1 - cor(mat, method="spearman")))
distM <- as.dist(refM) 
htree <- hclust(distM, method='ward.D2')
clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=200, distM=refM, verbose=0))
expect_identical(clusters, as.integer(obs))

obs <- quickCluster(mat, min.size=50)
clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=refM, verbose=0))
expect_identical(clusters, as.integer(obs))

mat <- matrix(rpois(10000, lambda=5), nrow=20)
subset.row <- 15:1 # With subsetting
refM <- sqrt(0.5*(1 - cor(mat[subset.row,], method="spearman")))
distM <- as.dist(refM) 
htree <- hclust(distM, method='ward.D2')
obs <- quickCluster(mat, min.size=50, subset.row=subset.row)
clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=refM, verbose=0))
expect_identical(clusters, as.integer(obs))

# Checking that we're executing the igraph methods correctly.

set.seed(300002)
test_that("quickCluster with igraph works with min.size settings", {
    mat <- matrix(rpois(10000, lambda=5), nrow=20)
    obs <- quickCluster(mat, min.size=0, method="igraph")
    ref <- quickCluster(mat, get.rank=TRUE)
    
    snn <- buildSNNGraph(ref) 
    library(igraph)
    out <- cluster_fast_greedy(snn)
    expect_identical(factor(out$membership), obs)
    
    min.size <- 100 # Checking that min.size merging works.
    expect_false(all(table(obs) >= min.size))
    obs2 <- quickCluster(mat, min.size=min.size, method="igraph")
    expect_true(all(table(obs2) >= min.size))

    combined <- paste0(obs, ".", obs2)
    expect_identical(length(unique(combined)), length(unique(obs))) # Confirm that they are nested.
})

# Other checks

expect_identical(length(quickCluster(mat, method="igraph", d=NA)), ncol(mat)) # Checking that the dimensions are correct for igraph.
suppressWarnings(expect_false(identical(quickCluster(mat), quickCluster(mat[subset.row,])))) # Checking that subsetting gets different results.
suppressWarnings(expect_identical(quickCluster(mat, subset.row=subset.row), quickCluster(mat[subset.row,]))) # Checking that subset.row works.

# Seeing how it interacts with the normalization method.

set.seed(300003)
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- seq_len(ngenes)/100
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0

out <- computeSumFactors(dummy, cluster=known.clusters)
expect_equal(out, count.sizes/mean(count.sizes)) # Even though there is a majority of DE, each pair of clusters is still okay.

out1 <- computeSumFactors(dummy, cluster=known.clusters, ref=1)
expect_equal(out, out1)
out2 <- computeSumFactors(dummy, cluster=known.clusters, ref=2)
expect_equal(out, out2)
out3 <- computeSumFactors(dummy, cluster=known.clusters, ref=3)
expect_equal(out, out3)

expect_error(computeSumFactors(dummy, cluster=known.clusters, ref=0), "'ref.clust' value not in 'clusters'")

# Checking out what happens with silly inputs.

expect_error(quickCluster(dummy[0,]), "rank variances of zero detected for a cell")
expect_error(quickCluster(dummy[,0]), "fewer cells than the minimum cluster size")

leftovers <- 100
expect_warning(forced <- quickCluster(dummy[,c(which(known.clusters==1), 
                                               which(known.clusters==2), 
                                               which(known.clusters==3)[seq_len(leftovers)])]), 
               sprintf("%i cells were not assigned to any cluster", leftovers))
expect_identical(as.character(tail(forced, leftovers)), rep("0", leftovers))

# Trying it out on a SCESet object.

set.seed(20002)
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- sample(seq_len(ngenes)/100)
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
emp.clusters <- quickCluster(X)
expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)

