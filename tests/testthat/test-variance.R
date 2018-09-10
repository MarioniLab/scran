# This tests the variance calculation functions in scran.
# require(scran); require(testthat); source("test-variance.R")

set.seed(20002)
true.p <- runif(100)
trended <- runif(100, 1, 2)
df <- 20
observed <- trended * qchisq(true.p, df=df, lower.tail=FALSE)/df

test_that("testVar's chi-squared test works as expected", {
    pvals <- testVar(observed, trended, df=df)
    expect_equal(pvals, true.p)
    pvals <- testVar(observed, trended, df=df, log.p=TRUE)
    expect_equal(pvals, log(true.p))
    
    design <- model.matrix(~factor(rep(c(1,2), each=11)))
    pvals <- testVar(observed, trended, design=design)
    expect_equal(pvals, true.p)
})

test_that("testVar's F-test works as expected", {
    nspikes <- ncells <- 100
    spike.means <- 2^runif(nspikes, -2, 8)
    spike.disp <- (100/spike.means + 0.5) * 1/rchisq(nspikes, df=10)
    spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)
    
    exprs <- log2(spike.data/(colSums(spike.data)/mean(colSums(spike.data)))+1)
    fit <- trendVar(exprs)
    pvals <- testVar(fit$var, fit$trend(fit$mean), df=ncells-1, second.df=fit$df2, test='f')
    
    rat <- (fit$var/fit$trend(fit$mean))
    df1 <- ncells - 1L
    ffit <- limma::fitFDistRobustly(rat[fit$var > 0 & fit$mean >= 0.1], df=df1) # filtering out zero-variance, low-abundance genes.
    expect_equal(pvals, pf(rat/ffit$scale, df1=df1, df2=ffit$df2, lower.tail=FALSE))

    # Testing log-transformation.
    lpvals <- testVar(fit$var, fit$trend(fit$mean), df=ncells-1, second.df=fit$df2, test='f', log.p=TRUE)
    expect_equal(log(pvals), lpvals)
    
    # Checking for error upon no df.2
    expect_error(testVar(fit$var, fit$trend(fit$mean), df=ncells-1, test='f'),
                 "second df from trendVar() must be specified for test='f'", fixed=TRUE)
})

test_that("testVar works with silly inputs", {
    expect_identical(testVar(0, 0, df=10), 1)
    expect_identical(testVar(0, 0, df=10, test="f", second.df=2), 1)
    expect_identical(testVar(0, 0, df=0), NA_real_)
    expect_identical(testVar(0, 0, df=0, test="f", second.df=2), NA_real_)

    expect_identical(testVar(numeric(0), trended, df=df), numeric(0))
    expect_identical(testVar(numeric(0), trended, df=df, test="f", second.df=2), numeric(0))

    expect_identical(testVar(observed, numeric(0), df=df), numeric(0))
    expect_identical(testVar(observed, numeric(0), df=df, test="f", second.df=2), numeric(0))

    expect_identical(testVar(observed, trended, df=numeric(0)), rep(NA_real_, length(observed)))
    expect_identical(testVar(observed, trended, df=numeric(0), test="f", second.df=2), rep(NA_real_, length(observed)))
})

####################################################################################################

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
                                    1, FUN=function(p) { pnorm(sum(qnorm(p) * DF)/sum(DF)) }))

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

    # Checking p-value behaviour with p-values at the boundary.
    dec$p.value[1] <- 0
    dec2$p.value[1] <- 0
    out <- combineVar(dec, dec2, method="z")
    expect_equal(out$p.value[1], 0)
    out <- combineVar(dec, dec2, method="berger")
    expect_equal(out$p.value[1], 0)
    out <- combineVar(dec, dec2, method="simes")
    expect_equal(out$p.value[1], 0)
    out <- combineVar(dec, dec2, method="fisher")
    expect_equal(out$p.value[1], 0)

    dec2$p.value[1] <- 1
    out <- combineVar(dec, dec2, method="z")
    expect_equal(out$p.value[1], 0.5)
    out <- combineVar(dec, dec2, method="berger")
    expect_equal(out$p.value[1], 1)
    out <- combineVar(dec, dec2, method="simes")
    expect_equal(out$p.value[1], 0)
    out <- combineVar(dec, dec2, method="fisher")
    expect_equal(out$p.value[1], 0)

    # Checking empty inputs.
    out <- combineVar(dec[0,], dec2[0,], dec3[0,])
    expect_equal(nrow(out), 0L)
    expect_identical(colnames(out), c("mean", "total", "bio", "tech", "p.value", "FDR"))
})
