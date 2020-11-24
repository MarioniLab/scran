# This tests the combineVar function.
# require(scran); require(testthat); source("test-combine-var.R")

set.seed(20003)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
X <- scuttle::logNormCounts(X)
dec <- modelGeneVar(X)

sub.d <- X[,seq_len(ncells/2)]
block <- sample(3, replace=TRUE, ncol(sub.d))
dec2 <- modelGeneVar(sub.d, block=block)

alt.d <- X[,ncells/2+1:50]
design <- model.matrix(~runif(ncol(alt.d)))    
dec3 <- modelGeneVar(alt.d, design=design)

test_that("combineVar works correctly", {
    # Checking averaging of stats.
    res <- combineVar(dec, dec2, dec3, method="stouffer")
    expect_equal(res$mean, rowMeans(cbind(dec$mean, dec2$mean, dec3$mean)))
    expect_equal(res$total, rowMeans(cbind(dec$total, dec2$total, dec3$total)))
    expect_equal(res$tech, rowMeans(cbind(dec$tech, dec2$tech, dec3$tech)))
    expect_equal(res$bio, rowMeans(cbind(dec$bio, dec2$bio, dec3$bio)))
    expect_identical(rownames(res), rownames(dec))

    # Checking proper calculation of combined p-values.
    pvalmat <- cbind(dec$p.value, dec2$p.value, dec3$p.value)
    expect_equal(res$p.value, apply(pvalmat, 1, FUN=function(p) { pnorm(sum(qnorm(p))/sqrt(3)) } ))

    res2 <- combineVar(dec, dec2, dec3, method="simes")
    expect_equal(res[,c("mean", "total", "tech", "bio")], res2[,c("mean", "total", "tech", "bio")])
    expect_equivalent(res2$p.value, apply(pvalmat, 1, FUN=function(p) { min(p.adjust(p, method="BH")) }))

    res3 <- combineVar(dec, dec2, dec3, method="berger")
    expect_equal(res[,c("mean", "total", "tech", "bio")], res3[,c("mean", "total", "tech", "bio")])
    expect_equal(res3$p.value, apply(pvalmat, 1, max))

    res4 <- combineVar(dec, dec2, dec3, method="fisher")
    expect_equal(res[,c("mean", "total", "tech", "bio")], res4[,c("mean", "total", "tech", "bio")])
    expect_equivalent(res4$p.value, pchisq(-2*rowSums(log(pvalmat)), df=6, lower.tail=FALSE))

    # Same results with a list of DF's.
    expect_identical(res, combineVar(list(dec, dec2, dec3), method="stouffer"))
    expect_identical(res, combineVar(list(dec, dec2), dec3, method="stouffer"))
    expect_identical(res, combineVar(dec, list(dec2), dec3, method="stouffer"))
})

test_that("combineVar works when weighting is turned on", {
    ref <- combineVar(dec, dec2, dec3, method="stouffer") 
    res <- combineVar(dec, dec2, dec3, method="stouffer", equiweight=FALSE)
    expect_equal(res, ref)

    N <- c(ncells, ncol(sub.d), ncol(alt.d))
    res <- combineVar(dec, dec2, dec3, method="stouffer", equiweight=FALSE, ncells=N)
    expect_equal(res$mean, drop(cbind(dec$mean, dec2$mean, dec3$mean) %*% N)/sum(N))
    expect_equal(res$bio, drop(cbind(dec$bio, dec2$bio, dec3$bio) %*% N)/sum(N))
    expect_equal(res$total, drop(cbind(dec$total, dec2$total, dec3$total) %*% N)/sum(N))
    expect_equal(res$tech, drop(cbind(dec$tech, dec2$tech, dec3$tech) %*% N)/sum(N))

    # Checking proper calculation of combined p-values.
    pvalmat <- cbind(dec$p.value, dec2$p.value, dec3$p.value)
    expect_equal(res$p.value, apply(pvalmat, 1, FUN=function(p) { pnorm(sum(N*qnorm(p))/sqrt(sum(N^2))) } ))

    # Other methods are unaffected.
    ref <- combineVar(dec, dec2, dec3, method="fisher")
    res <- combineVar(dec, dec2, dec3, method="fisher", equiweight=FALSE, ncells=N)
    expect_equivalent(res$p.value, ref$p.value)
})

test_that("combineVar behaves in edge cases", {
    # Just directly returns the input if only one DF is supplied.
    expect_equal(combineVar(dec), dec)
    expect_equal(combineVar(dec2), dec2)
    expect_equal(combineVar(dec3), dec3)

    # Checking failures:
    expect_error(res <- combineVar(dec, dec2[rev(rownames(dec)),]), "gene identities should be the same") 

    # Checking empty inputs.
    out <- combineVar(dec[0,], dec2[0,], dec3[0,])
    expect_equal(nrow(out), 0L)
    expect_identical(colnames(out), c("mean", "total", "tech", "bio", "p.value", "FDR", "per.block"))
})

#######################################
#######################################

dec <- modelGeneCV2(X)
dec2 <- modelGeneCV2(sub.d, block=block)
dec3 <- modelGeneCV2(alt.d)

geoRowMeans <- function(mat) {
    exp(rowMeans(log(mat)))
}

test_that("combineCV2 works correctly", {
    res <- combineCV2(dec, dec2, dec3)
    expect_equal(res$mean, geoRowMeans(cbind(dec$mean, dec2$mean, dec3$mean)))
    expect_equal(res$total, geoRowMeans(cbind(dec$total, dec2$total, dec3$total)))
    expect_equal(res$trend, geoRowMeans(cbind(dec$trend, dec2$trend, dec3$trend)))
    expect_equal(res$ratio, geoRowMeans(cbind(dec$ratio, dec2$ratio, dec3$ratio)))
    expect_equivalent(res$p.value, metapod::parallelFisher(list(dec$p.value, dec2$p.value, dec3$p.value))$p.value)

    # Same results with a list of DF's.
    expect_identical(res, combineCV2(list(dec, dec2, dec3)))
    expect_identical(res, combineCV2(list(dec, dec2), dec3))
    expect_identical(res, combineCV2(dec, list(dec2), dec3))
})

geoRowMeansW <- function(mat, w) {
    exp(drop(log(mat) %*% w / sum(w)))
}

test_that("combineCV2 works when weighting is turned on", {
    ref <- combineCV2(dec, dec2, dec3, method="stouffer") 
    res <- combineCV2(dec, dec2, dec3, method="stouffer", equiweight=FALSE)
    expect_equal(res, ref)

    N <- c(ncells, ncol(sub.d), ncol(alt.d))
    res <- combineCV2(dec, dec2, dec3, method="stouffer", equiweight=FALSE, ncells=N)
    expect_equal(res$mean, geoRowMeansW(cbind(dec$mean, dec2$mean, dec3$mean), N))
    expect_equal(res$ratio, geoRowMeansW(cbind(dec$ratio, dec2$ratio, dec3$ratio), N))
    expect_equal(res$total, geoRowMeansW(cbind(dec$total, dec2$total, dec3$total), N))
    expect_equal(res$trend, geoRowMeansW(cbind(dec$trend, dec2$trend, dec3$trend), N))
})

test_that("combineCV2 behaves in edge cases", {
    # Just directly returns the input if only one DF is supplied.
    expect_equal(combineCV2(dec), dec)
    expect_equal(combineCV2(dec2), dec2)
    expect_equal(combineCV2(dec3), dec3)

    # Checking failures:
    expect_error(res <- combineCV2(dec, dec2[rev(rownames(dec)),]), "gene identities should be the same") 

    # Checking empty inputs.
    out <- combineCV2(dec[0,], dec2[0,], dec3[0,])
    expect_equal(nrow(out), 0L)
    expect_identical(colnames(out), c("mean", "total", "trend", "ratio", "p.value", "FDR", "per.block"))
})
