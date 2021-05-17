# This provides some tests for scoreMarkers.
# library(testthat); library(scran); source("test-score-markers.R")

set.seed(10000)

split_by_cluster <- function(y, cluster) {
    by.clust <- split(seq_along(cluster), cluster)
    lapply(by.clust, function(i) y[,i,drop=FALSE])
}

test_that("Cohen calculations are collated correctly", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- sample(letters[1:5], ncol(y), replace=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 
    mats <- split_by_cluster(y, cluster)

    for (i in names(mats)) {
        current <- out[[i]]$full.logFC.cohen
        lmat <- mats[[i]]
        left <- rowMeans(lmat)
        lvar <- rowVars(lmat)

        for (j in setdiff(names(mats), i)) {
            rmat <- mats[[j]]
            right <- rowMeans(rmat)
            rvar <- rowVars(rmat)

            expect_equal(current[,j], (left - right)/sqrt( (lvar + rvar)/2 ))
            expect_equal(mcols(current)[j,], ncol(lmat) * ncol(rmat))

            flip <- out[[j]]$full.logFC.cohen[,i]
            expect_equal(flip, -current[,j])
        }
    }
})

test_that("Cohen calculations handle zero variances gracefully", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- rep(letters[1:2], length.out=ncol(y))
    ref <- scoreMarkers(y, cluster, full.stats=TRUE) 

    y[c(1, 10),cluster=='b'] <- 0
    y[c(1, 10),cluster=='a'] <- 1
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 

    for (i in names(ref)) {
        curref <- ref[[i]]$full.logFC.cohen
        curout <- out[[i]]$full.logFC.cohen
        curref[c(1,10),] <- Inf * (-1)^(i=="b")
        expect_identical(curref, curout)
    }

    # Special case if the log-fold change is zero.
    y[c(1, 10),cluster=='a'] <- 0
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 

    for (i in names(ref)) {
        curref <- ref[[i]]$full.logFC.cohen
        curout <- out[[i]]$full.logFC.cohen
        curref[c(1,10),] <- 0
        expect_identical(curref, curout)
    }
})

test_that("AUC effect size calculations are collated correctly", {
    y <- matrix(rpois(1000, 0.5), ncol=100)
    cluster <- sample(letters[1:5], ncol(y), replace=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 
    mats <- split_by_cluster(y, cluster)

    AUC_calculator <- function(left, right) { 
        U <- vapply(seq_len(nrow(left)), function(i) wilcox.test(left[i,], right[i,], exact=FALSE)$statistic, 0)
        U / ncol(left) / ncol(right)
    }

    for (i in names(mats)) {
        current <- out[[i]]$full.AUC
        lmat <- mats[[i]]

        for (j in setdiff(names(mats), i)) {
            rmat <- mats[[j]]
            ref <- AUC_calculator(lmat, rmat)
            expect_equal(current[,j], ref)

            flip <- out[[j]]$full.AUC[,i]
            expect_equal(flip, 1-current[,j])
        }
    }
})

test_that("AUC calculations handle ties and zeroes properly", {
    y <- matrix(0, ncol=100, nrow=5)
    cluster <- sample(letters[1:5], ncol(y), replace=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 

    for (i in names(out)) {
        expect_true(all(abs(as.matrix(out[[i]]$full.AUC) - 0.5) < 1e-8))
    }

    # Same behavior for non-zero tied values.
    y <- matrix(1, ncol=100, nrow=5)
    out2 <- scoreMarkers(y, cluster, full.stats=TRUE)

    for (i in names(out)) {
        expect_identical(out[[i]]$full.AUC, out2[[i]]$full.AUC)
    }
})

test_that("AUC calculations handles negative values properly", {
    y <- matrix(rpois(1000, 0.5), ncol=100)
    cluster <- sample(letters[1:5], ncol(y), replace=TRUE)

    ref <- scoreMarkers(y, cluster, full.stats=TRUE) 
    out <- scoreMarkers(-y, cluster, full.stats=TRUE) 

    for (i in names(out)) {
        expect_equal(as.matrix(out[[i]]$full.AUC), 1-as.matrix(ref[[i]]$full.AUC))
    }
})

test_that("logFC detected calculations are done properly", {
    y <- matrix(rbinom(1000, 1, 0.2), ncol=100)
    cluster <- sample(letters[1:5], ncol(y), replace=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 

    mats <- split_by_cluster(y, cluster)

    for (i in names(mats)) {
        current <- out[[i]]$full.logFC.detected
        ref.left <- rowMeans(mats[[i]] > 0)

        for (j in setdiff(names(mats), i)) {
            ref.right <- rowMeans(mats[[j]] > 0)
            lfc <- current[,j]

            # Shrinkage works correctly.
            pos <- lfc > 1e-8
            expect_identical(pos, ref.left > ref.right + 1e-8)
            expect_true(all(lfc[pos] < log2(ref.left[pos] / ref.right[pos])))

            neg <- lfc < -1e-8
            expect_identical(neg, ref.left < ref.right - 1e-8)
            expect_true(all(lfc[neg] > log2(ref.left[neg] / ref.right[neg])))
        }
    }
})

test_that("logFC detected are properly set to zero", {
    y0 <- matrix(rbinom(1000, 1, 0.2), ncol=20)
    y <- cbind(y0, y0)
    cluster <- rep(1:2, each=ncol(y0))
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 

    for (i in names(out)) {
        current <- out[[i]]$full.logFC.detected
        expect_true(all(as.matrix(current)==0))
    }
})

test_that("logFC detected avoids undefined log-fold changes", {
    y0 <- matrix(rbinom(1000, 1, 0.2), ncol=20)
    y <- cbind(y0, matrix(0, nrow(y0), ncol(y0)))
    cluster <- rep(1:2, each=ncol(y0))
    out <- scoreMarkers(y, cluster, full.stats=TRUE) 

    for (i in names(out)) {
        current <- out[[i]]$full.logFC.detected
        expect_true(all(is.finite(as.matrix(current))))
        if (i=="1") {
            expect_true(all(as.matrix(current) > 0))
        } else {
            expect_true(all(as.matrix(current) < 0))
        }
    }
})

test_that("effect summary statistics are correctly computed", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- rep(letters[5:10], length.out=ncol(y))
    ref <- scoreMarkers(y, cluster, full.stats=TRUE) 

    for (i in names(ref)) {
        current <- ref[[i]]

        for (n in c("logFC.cohen", "AUC", "logFC.detected")) {
            mat <- as.matrix(current[,paste0("full.", n)])
            expect_equal(current[,paste0("mean.", n)], rowMeans(mat))
            expect_equal(current[,paste0("max.", n)], rowMaxs(mat))
            expect_equal(current[,paste0("max.", n)], rowMaxs(mat))
            expect_equal(current[,paste0("median.", n)], rowMedians(mat))
            expect_equal(current[,paste0("rank.", n)], computeMinRank(mat))
        }
    }
})

test_that("non-effect summary statistics are correctly computed", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- rep(letters[5:10], length.out=ncol(y))
    ref <- scoreMarkers(y, cluster, full.stats=TRUE) 

    mats <- split_by_cluster(y, cluster)
    means <- do.call(cbind, lapply(mats, rowMeans))
    detect <- do.call(cbind, lapply(mats, function(x) rowMeans(x > 0)))

    for (i in names(ref)) {
        current <- ref[[i]]

        expect_equal(current$self.average, means[,i])
        expect_equal(current$other.average, rowMeans(means[,setdiff(colnames(means), i)]))

        expect_equal(current$self.detected, detect[,i])
        expect_equal(current$other.detected, rowMeans(detect[,setdiff(colnames(detect), i)]))
    }
})

test_that("blocking is performed correctly", {
    y <- matrix(rpois(1000, 0.5), ncol=100)
    cluster <- sample(letters[1:5], ncol(y), replace=TRUE)
    block <- sample(1:3, ncol(y), replace=TRUE)

    b1 <- scoreMarkers(y[,block==1], cluster[block==1], full.stats=TRUE)
    b2 <- scoreMarkers(y[,block==2], cluster[block==2], full.stats=TRUE)
    b3 <- scoreMarkers(y[,block==3], cluster[block==3], full.stats=TRUE)

    out <- scoreMarkers(y, cluster, full.stats=TRUE, block=block) 

    for (i in letters[1:5]) {
        # Effect size is computed correctly.
        current1 <- b1[[i]]$full.logFC.cohen
        current2 <- b2[[i]]$full.logFC.cohen
        current3 <- b3[[i]]$full.logFC.cohen

        for (j in setdiff(letters[1:5], i)) {
            vmat <- cbind(current1[,j], current2[,j], current3[,j])

            w1 <- mcols(current1)[j,]
            w2 <- mcols(current2)[j,]
            w3 <- mcols(current3)[j,]
            wmat <- t(t(!is.na(vmat)) * c(w1, w2, w3))

            expect_equal(out[[i]]$full.logFC.cohen[,j], rowSums(vmat * wmat)/rowSums(wmat))
        }
    }

    # Summary stats make sense.
    for (i in names(out)) {
        current <- out[[i]]
        summed <- current$self.average + current$other.average * 4L
        if (i==names(out)[1]) {
            ref <- summed
            first <- current$self.average
        } else {
            expect_equal(summed, ref)
            expect_false(isTRUE(all.equal(first, current$self.average)))
        }
    }
})

test_that("blocking is performed correctly (group-specific batches)", {
    y <- matrix(rpois(1000, 0.5), ncol=100)
    cluster <- sample(10:14, ncol(y), replace=TRUE)
    block <- ifelse(cluster %% 2 == 1, "odd", "even")

    bo <- scoreMarkers(y[,block=="odd"], cluster[block=="odd"], full.stats=TRUE)
    be <- scoreMarkers(y[,block=="even"], cluster[block=="even"], full.stats=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE, block=block) 

    onames <- as.character(10:14)
    for (i in onames) {
        currentO <- bo[[i]]$full.logFC.detected
        currentE <- be[[i]]$full.logFC.detected

        for (j in setdiff(onames, i)) {
            if ((as.integer(j) %% 2) != (as.integer(i) %% 2)) {
                expect_equal(out[[i]]$full.logFC.detected[,j], rep(NA_real_, nrow(out[[i]])))
            } else {
                if (as.integer(j) %% 2) {
                    vmat <- currentO[,j]
                } else {
                    vmat <- currentE[,j]
                }
                expect_equal(out[[i]]$full.logFC.detected[,j], vmat)
            }
        }
    }

    # Summary stats are valid.
    for (i in names(out)) {
        current <- out[[i]][,1:4]
        expect_true(all(is.finite(as.matrix(current))))
    }
})

test_that("blocking is performed correctly (batch-specific groups)", {
    y <- matrix(rpois(1000, 0.5), ncol=100)
    cluster <- sample(10:14, ncol(y), replace=TRUE)
    block <- sample(2, ncol(y), replace=TRUE)
    block[cluster==10] <- 1
    block[cluster==14] <- 2

    b1 <- scoreMarkers(y[,block==1], cluster[block==1], full.stats=TRUE)
    b2 <- scoreMarkers(y[,block==2], cluster[block==2], full.stats=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE, block=block) 

    is.shared <- !cluster %in% c(10, 14)
    shared <- scoreMarkers(y[,is.shared], cluster[is.shared], full.stats=TRUE, block=block[is.shared])

    onames <- as.character(10:14)
    for (i in onames) {
        current1 <- b1[[i]]$full.AUC
        current2 <- b2[[i]]$full.AUC

        for (j in setdiff(onames, i)) {
            if (i=="10") {
                if (j!="14") {
                    ref <- current1[,j]
                } else {
                    ref <- rep(NA_real_, nrow(y))
                }
            } else if (i=="14") {
                if (j!="10") {
                    ref <- current2[,j]
                } else {
                    ref <- rep(NA_real_, nrow(y))
                }
            } else {
                if (j=="10") {
                    ref <- current1[,j]
                } else if (j=="14") {
                    ref <- current2[,j]
                } else {
                    ref <- shared[[i]]$full.AUC[,j]
                }
            }
            expect_equal(ref, out[[i]]$full.AUC[,j])
        }
    }

    # Summary stats are valid.
    for (i in names(out)) {
        current <- out[[i]][,1:4]
        expect_true(all(is.finite(as.matrix(current))))
    }
})

set.seed(100002)
test_that("group ordering is correct", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- sample(1:5, ncol(y), replace=TRUE)
    cluster[1:5] <- c(1,5,2,4,3) # guaranteeing that the first few elements are out of order.
    ref <- scoreMarkers(y, cluster, full.stats=TRUE) 

    n <- as.character(1:5)
    expect_identical(names(ref), n)
    for (i in n) {
        expect_identical(colnames(ref[[i]]$full.AUC), setdiff(n, i))
    }

    super_unname <- function(x) {
        x <- as.matrix(x)
        colnames(x) <- NULL
        x
    }

    out <- scoreMarkers(y, LETTERS[cluster], full.stats=TRUE) 
    n <- LETTERS[1:5]
    expect_identical(names(out), n)
    for (i in seq_along(n)) {
        expect_identical(colnames(out[[i]]$full.AUC), setdiff(n, n[i]))
        expect_identical(super_unname(ref[[i]]$full.AUC), super_unname(out[[i]]$full.AUC))
    }

    # Handles factors with unordered factors.
    n <- factor(LETTERS[1:5], LETTERS[5:1])
    out2 <- scoreMarkers(y, n[cluster], full.stats=TRUE) 
    expect_identical(names(out2), levels(n))
    for (i in LETTERS[1:5]) {
        past <- out[[i]]$full.AUC
        curr <- out2[[i]]$full.AUC
        expect_equal(past, curr[,rev(colnames(curr))])
    }

    # Handles factors with unused levels.
    n <- factor(LETTERS[2:6], LETTERS[1:10])
    ref <- scoreMarkers(y, n[cluster], full.stats=TRUE) 
    out <- scoreMarkers(y, as.character(n)[cluster], full.stats=TRUE) 
    expect_identical(ref, out)
})

set.seed(100002)
test_that("lfc handling is correct", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- sample(c(1,"A",3,"C",5), ncol(y), replace=TRUE)
    out <- scoreMarkers(y, cluster, full.stats=TRUE, lfc=1.5)
    
    for (x in c("1", "A", "3", "C", "5")) {
        y2 <- y
        y2[,cluster==x] <- y2[,cluster==x] - 1.5
        ref <- scoreMarkers(y2, cluster, full.stats=TRUE)

        expect_equal(out[[x]]$full.logFC.cohen, ref[[x]]$full.logFC.cohen)
        expect_equal(out[[x]]$mean.logFC.cohen, ref[[x]]$mean.logFC.cohen)
        expect_equal(out[[x]]$mean.AUC, ref[[x]]$mean.AUC)
        expect_equal(out[[x]]$full.AUC, ref[[x]]$full.AUC)
    }

    # Detected log-fold change is unaffected.
    ref2 <- scoreMarkers(y, cluster, full.stats=TRUE)

    for (x in c("1", "A", "3", "C", "5")) {
        expect_equal(out[[x]]$mean.logFC.detected, ref2[[x]]$mean.logFC.detected)
        expect_equal(out[[x]]$full.logFC.detected, ref2[[x]]$full.logFC.detected)
        expect_equal(out[[x]][,1:4], ref2[[x]][,1:4])
    }
})

set.seed(100003)
test_that("row.data handling is correct", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- sample(c(1,"A",3,"C",5), ncol(y), replace=TRUE)
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    rd <- DataFrame(row.names=rownames(y), stuff=runif(nrow(y)))

    ref <- scoreMarkers(y, cluster)
    out <- scoreMarkers(y, cluster, row.data=rd)
    for (i in seq_along(ref)) {
        expect_identical(out[[i]][,1], rd$stuff)
        expect_identical(out[[i]]$stuff, rd$stuff)
        expect_identical(ref[[i]]$self.average, out[[i]]$self.average)
    }
})

set.seed(100003)
test_that("subset handling is correct", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- sample(c(1,"A",3,"C",5), ncol(y), replace=TRUE)

    sub <- sample(nrow(y), 5)
    ref <- scoreMarkers(y[sub,], cluster)
    out <- scoreMarkers(y, cluster, subset.row=sub)
    expect_identical(ref, out)
    
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    csub <- rownames(y)[sub]
    ref <- scoreMarkers(y[csub,], cluster)
    out <- scoreMarkers(y, cluster, subset.row=csub)
    expect_identical(ref, out)

    # Behaves correctly with row.data
    rd <- DataFrame(row.names=rownames(y), stuff=runif(nrow(y)))
    ref <- scoreMarkers(y[csub,], cluster, row.data=rd[csub,,drop=FALSE])
    out <- scoreMarkers(y, cluster, subset.row=csub, row.data=rd)
    expect_identical(ref, out)
})

set.seed(100004)
test_that("pairings handling is correct", {
    y <- matrix(rnorm(1000), ncol=100)
    cluster <- sample(c(1,"A",3,"C",5), ncol(y), replace=TRUE)

    # Equivalent to subsetting the groups.
    keep <- cluster %in% as.character(1:5)
    ref <- scoreMarkers(y[,keep], cluster[keep])
    out <- scoreMarkers(y, cluster, pairings=as.character(1:5))

    expect_identical(ref, out)

    # A more directed approach.
    mat <- rbind(c("A", "C"), c(1,5))
    out <- scoreMarkers(y, cluster, pairings=mat)
    expect_identical(sort(names(out)), sort(mat[,1]))

    out1 <- scoreMarkers(y[,cluster %in% mat[1,]], cluster[cluster %in% mat[1,]])
    expect_identical(out[["A"]], out1[["A"]])

    out2 <- scoreMarkers(y[,cluster %in% mat[2,]], cluster[cluster %in% mat[2,]])
    expect_identical(out[["1"]], out2[["1"]])
})

