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
X <- normalize(X)
d <- exprs(X)
out <- trendVar(d)
dec <- decomposeVar(d, out)

sub.d <- d[,seq_len(ncells/2)]
block <- sample(3, replace=TRUE, ncol(sub.d))
out2 <- trendVar(sub.d, block=block)
dec2 <- decomposeVar(sub.d, out2)

alt.d <- d[,ncells/2+1:50]
design <- model.matrix(~runif(ncol(alt.d)))    
out3 <- trendVar(alt.d, design=design)
dec3 <- decomposeVar(alt.d, out3)

test_that("combineVar works correctly", {
    # Checking averaging of stats.
    N <- c(ncells, ncol(sub.d), ncol(alt.d))
    DF <- c(ncells - 1L, ncol(sub.d) - 3L, ncol(alt.d) - 2L)
    res <- combineVar(dec, dec2, dec3, method="z")
    expect_equal(res$mean, drop(cbind(dec$mean, dec2$mean, dec3$mean) %*% N / sum(N)))
    expect_equal(res$total, drop(cbind(dec$total, dec2$total, dec3$total) %*% DF) / sum(DF))
    expect_equal(res$tech, drop(cbind(dec$tech, dec2$tech, dec3$tech) %*% DF) / sum(DF))
    expect_equal(res$bio, drop(cbind(dec$bio, dec2$bio, dec3$bio) %*% DF) / sum(DF))

    # Checking proper calculation of combined p-values.
    expect_equal(res$p.value, apply(cbind(dec$p.value, dec2$p.value, dec3$p.value),
                                    1, FUN=function(p) { pnorm(sum(qnorm(p) * DF)/sqrt(sum(DF^2))) }))

    res2 <- combineVar(dec, dec2, dec3, method="simes")
    expect_equal(res[,c("mean", "total", "tech", "bio")], res2[,c("mean", "total", "tech", "bio")])
    expect_equal(res2$p.value, apply(cbind(dec$p.value, dec2$p.value, dec3$p.value),
                                     1, FUN=function(p) { min(p.adjust(p, method="BH")) }))

    res3 <- combineVar(dec, dec2, dec3, method="berger")
    expect_equal(res[,c("mean", "total", "tech", "bio")], res3[,c("mean", "total", "tech", "bio")])
    expect_equal(res3$p.value, apply(cbind(dec$p.value, dec2$p.value, dec3$p.value), 1, max))

    res4 <- combineVar(dec, dec2, dec3, method="fisher")
    expect_equal(res[,c("mean", "total", "tech", "bio")], res4[,c("mean", "total", "tech", "bio")])
    expect_equal(res4$p.value, pchisq(-2*rowSums(log(cbind(dec$p.value, dec2$p.value, dec3$p.value))),
                                      df=6, lower.tail=FALSE))
})

test_that("combineVar responds to settings", {
    res <- combineVar(dec, dec2, dec3)

    # Same results upon subsetting.
    reres <- combineVar(dec[1:10,], dec2[1:10,], dec3[1:10,])
    rescheck <- res[1:10,]
    rescheck$FDR <- p.adjust(rescheck$p.value, method="BH")
    expect_equal(reres, rescheck)

    # Just directly returns the input if only one DF is supplied.
    expect_equal(combineVar(dec), dec)
    expect_equal(combineVar(dec2), dec2)
    expect_equal(combineVar(dec3), dec3)

    # Checking that it performs correctly without weighting.
    res.unw <- combineVar(dec, dec2, dec3, weighted=FALSE)
    expect_equal(metadata(res.unw), metadata(res))

    dec.x <- dec
    dec2.x <- dec2
    dec3.x <- dec3
    metadata(dec.x) <- metadata(dec2.x) <- metadata(dec3.x) <- list(num.cells=1, resid.df=1)
    res.unw.ref <- combineVar(dec.x, dec2.x, dec3.x, weighted=TRUE)
    metadata(res.unw.ref) <- metadata(res.unw)
    expect_equal(res.unw, res.unw.ref)
})

test_that("combineVar handles edge cases properly", {
    # Checking failures:
    dec3.x <- dec3
    metadata(dec3.x) <- list()
    expect_error(res <- combineVar(dec, dec3.x), "inputs should come from decomposeVar()", fixed=TRUE)
    expect_error(res <- combineVar(dec, dec2[rev(rownames(dec)),]), "gene identities should be the same") # when you switch up the order.

    # Checking empty inputs.
    out <- combineVar(dec[0,], dec2[0,], dec3[0,])
    expect_equal(nrow(out), 0L)
    expect_identical(colnames(out), c("mean", "total", "bio", "tech", "p.value", "FDR"))
})
