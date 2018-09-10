# This tests the variance calculation functions in scran.
# require(scran); require(testthat); source("test-variance.R")

set.seed(20001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
isSpike(X, "MySpikes") <- rbinom(ngenes, 1, 0.7)==0
sizeFactors(X) <- colSums(dummy)
sizeFactors(X, "MySpikes") <- runif(ncells)
X <- normalize(X)
fit <- trendVar(X)

test_that("Variance decomposition is working correctly", {
    # Testing that the metrics are all calculated correctly.
    out.all <- decomposeVar(X, fit, get.spikes=TRUE)
    ref.mean <- rowMeans(exprs(X))
    expect_equivalent(out.all$mean, ref.mean)
    ref.var <- apply(exprs(X), 1, var)
    expect_equivalent(ref.var, out.all$total)
    expect_equivalent(out.all$tech, fit$trend(ref.mean))
    expect_equivalent(out.all$bio, out.all$total-out.all$tech)

    # Checking that stats are correctly stored.
    expect_equal(metadata(out.all)$num.cells, ncells)
    expect_equal(metadata(out.all)$resid.df, ncells-1L)

    # Checking that the p-value calculation is correct.
    ref.p <- testVar(out.all$total, out.all$tech, df=ncells-1)
    expect_equivalent(ref.p, out.all$p.value)
    expect_equivalent(p.adjust(ref.p, method="BH"), out.all$FDR)

    # Testing that spike-in selection is working.
    out <- decomposeVar(X, fit, get.spikes=FALSE)
    ref <- decomposeVar(X[!isSpike(X),], fit, get.spikes=TRUE)
    expect_identical(out, ref)

    ref2 <- decomposeVar(X, fit, get.spikes=NA)
    out.all$p.value[isSpike(X)] <- NA
    out.all$FDR <- p.adjust(out.all$p.value, method="BH")
    expect_identical(ref2, out.all)

    expect_equal(metadata(out)$num.cells, ncells)
    expect_equal(metadata(out)$resid.df, ncells-1L)
})
   
test_that("decomposeVar behaves correctly with subsetting", {
    shuffled <- c(500:1, 501:1000)
    out.ref <- decomposeVar(X[shuffled,], fit)
    out2 <- decomposeVar(X, fit, subset.row=shuffled)
    expect_identical(out.ref, out2)  

    # Proper interaction with get.spikes.
    out.ref2 <- decomposeVar(exprs(X)[shuffled,], fit)
    was.spike <- which(isSpike(X)[shuffled])
    out.ref2$p.value[was.spike] <- NA
    out.ref2$FDR <- p.adjust(out.ref$p.value, method="BH")
    expect_identical(out.ref2, out2)

    out3 <- decomposeVar(X, fit, subset.row=shuffled, get.spikes=FALSE)
    out.ref3 <- decomposeVar(exprs(X)[setdiff(shuffled, which(isSpike(X))),], fit)
    expect_identical(out.ref3, out3)

    # Checks what happens when interaction with get.spikes is disabled.
    out.ref <- decomposeVar(X[shuffled,], fit, get.spikes=TRUE)
    out2 <- decomposeVar(X, fit, subset.row=shuffled, get.spikes=TRUE)
    expect_identical(out.ref, out2) 
    out.ref2 <- decomposeVar(exprs(X)[shuffled,], fit)
    expect_identical(out.ref2, out2)
})

test_that("decomposeVar checks size factor centering", {
    subX <- X[,1:10] # Checking that it raises a warning upon subsetting (where the size factors are no longer centered).
    expect_warning(decomposeVar(subX, fit), "size factors not centred")
    expect_warning(decomposeVar(normalize(subX), fit), NA)
})

test_that("decomposeVar works with all genes", {
    # Using all genes for trend fitting.
    all.fit <- trendVar(X, use.spikes=NA)
    all.dec <- decomposeVar(X, all.fit, get.spikes=TRUE)
    expect_equal(all.fit$mean, setNames(all.dec$mean, rownames(all.dec)))
    expect_equal(all.fit$var, setNames(all.dec$total, rownames(all.dec)))
    expect_equal(all.fit$trend(all.fit$mean), setNames(all.dec$tech, rownames(all.dec)))

    # Using only the elements in the fit object, and checking we get the same result.
    all.dec2 <- decomposeVar(fit=all.fit)
    expect_equal(all.dec, all.dec2)

    # Checking that it ignores any other value of 'block' or 'design' or 'subset.row'.
    blocking <- rep(1:3, length.out=ncol(X))
    expect_equal(all.dec2, decomposeVar(fit=all.fit, block=blocking))
    expect_equal(all.dec2, decomposeVar(fit=all.fit, design=model.matrix(~factor(blocking))))
    expect_equal(all.dec2, decomposeVar(fit=all.fit, subset.row=1:10))

    # Specifying 'fit' also works with blocking.
    all.fit.b <- trendVar(X, use.spikes=NA, block=blocking)
    all.dec.b <- decomposeVar(X, all.fit.b, get.spikes=TRUE)
    all.dec.b2 <- decomposeVar(fit=all.fit.b)
    expect_equal(all.dec.b, all.dec.b2)
})

test_that("decomposeVar works with design matrices", {
    # Testing with a trivial design matrix.
    design0 <- matrix(1, ncol(X), 1)
    fit <- trendVar(X)
    out <- decomposeVar(X, fit, get.spikes=FALSE)
    fit2 <- trendVar(X, design=design0)
    out2 <- decomposeVar(X, fit2, get.spikes=FALSE) # defaults to all-ones.
    expect_equal(out, out2)

    # Testing with a modified design matrix.
    design <- model.matrix(~factor(rep(c(1,2), each=100)))
    fit3 <- trendVar(X, design=design)
    out3 <- decomposeVar(X, fit3, get.spikes=FALSE)
    expect_equal(out$mean, out3$mean)

    refit <- lm.fit(y=t(exprs(X)), x=design)
    effects <- refit$effects[-seq_len(ncol(design)),]
    test.var <- colMeans(effects^2)
    
    expect_equivalent(out3$total, test.var[!isSpike(X)])
    expect_equivalent(out3$tech, fit3$trend(out$mean))
    expect_equivalent(out3$bio, out3$total-out3$tech)
    
    ref.p <- testVar(out3$total, out3$tech, df=nrow(design) - ncol(design))
    expect_equivalent(ref.p, out3$p.value)
    expect_equivalent(p.adjust(ref.p, method="BH"), out3$FDR)

    # Checking what happens when I use design in decomposeVar only.
    out2a <- decomposeVar(X, fit, design=design, get.spikes=FALSE)
    expect_equivalent(out2a$mean, out$mean)
    expect_equivalent(out2a$total, test.var[!isSpike(X)])
    expect_equivalent(out2a$tech, fit$trend(out$mean))
    expect_equivalent(out2a$bio, out2a$total-out2a$tech)
})

test_that("decomposeVar works with blocking", {
    # Testing it out.
    block0 <- rep(1, ncol(X))
    fit <- trendVar(X)
    out <- decomposeVar(X, fit)
    fit2 <- trendVar(X, block=block0)
    out2 <- decomposeVar(X, fit2)
    expect_equal(out, out2)
    
    # Trying it out with actual blocking, and manually checking the outputs.
    block <- sample(3, ncol(X), replace=TRUE)
    fit3 <- trendVar(X, block=block)
    out3 <- decomposeVar(X, fit3)

    var.mat <- cbind(apply(logcounts(X)[,block==1], 1, var),
                     apply(logcounts(X)[,block==2], 1, var),
                     apply(logcounts(X)[,block==3], 1, var))
    resid.df <- tabulate(block) - 1L
    total.var <- as.numeric(var.mat %*% resid.df/sum(resid.df))
    expect_equal(total.var, out3$total)

    total.mean <- unname(rowMeans(logcounts(X)))
    expect_equal(total.mean, out3$mean)

    mean.mat <- cbind(rowMeans(logcounts(X)[,block==1]),
                      rowMeans(logcounts(X)[,block==2]),
                      rowMeans(logcounts(X)[,block==3]))
    tech.var <- fit3$trend(as.vector(mean.mat))
    dim(tech.var) <- dim(mean.mat)
    total.tech <- as.numeric(tech.var %*% resid.df/sum(resid.df))
    expect_equal(total.tech, out3$tech)
    expect_equal(total.var - total.tech, out3$bio)
    
    pval <- testVar(as.numeric(var.mat), null=as.numeric(tech.var), 
                    df=rep(resid.df, each=nrow(X)))
    dim(pval) <- dim(mean.mat)
    pval[isSpike(X),] <- NA
    expect_equal(out3$p.value, pchisq(rowSums(log(pval))*-2, df=2*ncol(mean.mat), lower.tail=FALSE))

    # Checking out what happens when I use a normal trend but use block= in decomposeVar only.
    out2a <- decomposeVar(X, fit, block=block)
    expect_equal(total.var, out2a$total)
    expect_equal(total.mean, out2a$mean)

    tech.var <- fit$trend(as.vector(mean.mat))
    dim(tech.var) <- dim(mean.mat)
    total.tech <- as.numeric(tech.var %*% resid.df/sum(resid.df))
    expect_equal(total.tech, out2a$tech)
    expect_equal(total.var - total.tech, out2a$bio)
    
    # Checking what happens with a no-residual level.
    block[1] <- 0
    fit4a <- trendVar(X, block=block)
    out4a <- decomposeVar(X, fit4a)
    fit4b <- trendVar(X[,-1], block=block[-1])
    out4b <- decomposeVar(X[,-1], fit4b)
    expect_equal(out4a$total, out4b$total)
    expect_equal(out4a$tech, out4b$tech)
    expect_equal(out4a$bio, out4b$bio)
    expect_equal(out4a$mean, total.mean)
})

####################################################################################################

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
