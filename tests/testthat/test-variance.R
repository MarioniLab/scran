# This tests the variance calculation functions in scran.
# require(scran); require(testthat); source("test-variance.R")

set.seed(20000)
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

test_that("trendVar works on a basic scenario", {
    expect_equal(out$mean, rowMeans(d))
    expect_equal(out$var, apply(d, 1, var))
    expect_identical(out$design, NULL)
    expect_identical(out$design, NULL)
   
    # Hard to test it without copying the code, so I'll just check that the bounds are right.
    expect_is(out$trend, "function") 
    mx <- max(out$mean)
    mn <- min(out$mean)
    expect_equal(out$trend(mx), out$trend(mx+1))
    expect_equal(out$trend(0), 0)
    expect_equal(out$trend(1:10), sapply(1:10, out$trend)) # Checking we get consistent results with returned function.
    expect_equal(out$trend(100:1/20), sapply(100:1/20, out$trend))  # More checking, reversed order.
})

test_that("trendVar subsets properly (user-supplied and automatic)", {
    # Checking that genes with no variance don't break the results, but still get reported in the output.
    dz <- rbind(1, d, 0)
    outz <- trendVar(dz)
    expect_equal(outz$mean, c(1, out$mean, 0))
    expect_equal(outz$var, c(0, out$var, 0))
    expect_equal(out$trend(outz$mean), outz$trend(outz$mean))

    # Checking that it properly ignores the low-abundance genes.
    filt <- trendVar(d, min.mean=1)
    expect_equal(filt$mean, out$mean)
    expect_equal(filt$var, out$var)
    expect_false(all(rowMeans(d)>=1))
    ref <- trendVar(d[rowMeans(d)>=1,])
    expect_equal(filt$trend(0:100/10), ref$trend(0:100/10))

    # Checking that manual subset.row works. 
    shuffled <- c(500:110)
    out.ref <- trendVar(d[shuffled,])
    out2 <- trendVar(d, subset.row=shuffled)
    expect_equal(out.ref, out2) 
})

test_that("trendVar works correctly on SingleCellExperiment objects", {
    suppressWarnings(expect_error(trendVar(X, parametric=TRUE), "no spike-in transcripts present for 'use.spikes=TRUE'"))
    suppressWarnings(expect_error(trendVar(X, parametric=FALSE), "no spike-in transcripts present for 'use.spikes=TRUE'"))

    cntrl_data <- list(All=!logical(ngenes), Some=rbinom(ngenes, 1, 0.5)==0, None=logical(ngenes))
    test.x <- runif(1000, 0, max(out$mean))

    # Checking that spike-in referencing works.
    isSpike(X, "All") <- cntrl_data$All
    expect_identical(isSpike(X), cntrl_data$All) # Just checking here...
    out2 <- trendVar(X)
    expect_equal(out$mean, out2$mean)
    expect_equal(out$var, out2$var)
    expect_equal(out$trend(test.x), out2$trend(test.x))
    expect_equal(out$design, out2$design)
    
    # Checking that spike-in referencing fails, and use.spikes=FALSE works.
    isSpike(X, "All") <- NULL
    isSpike(X, "None") <- cntrl_data$None
    expect_identical(isSpike(X), cntrl_data$None)
    expect_error(trendVar(X), "no spike-in transcripts present for 'use.spikes=TRUE'")
    out3 <- trendVar(X, use.spikes=FALSE)
    expect_equal(out3$mean, out2$mean)
    expect_equal(out3$var, out2$var)
    expect_equal(out3$trend(test.x), out2$trend(test.x))
    expect_equal(out3$design, out2$design)
    
    # Checking that subsetting by spike-ins works.
    isSpike(X, "None") <- NULL
    isSpike(X, "Some") <- cntrl_data$Some
    expect_identical(isSpike(X), cntrl_data$Some)
    out3a <- trendVar(X)
    expect_equal(out3$mean[cntrl_data$Some], out3a$mean)
    expect_equal(out3$var[cntrl_data$Some], out3a$var)
    out3b <- trendVar(X, use.spikes=NA)
    expect_equal(out3$mean, out3b$mean)
    expect_equal(out3$var, out3b$var)
    expect_equal(out3$trend(test.x), out3b$trend(test.x))
    expect_equal(out3$design, out3b$design)
   
    # Checking for proper interaction between use.spike and subset.row.
    out3c <- trendVar(X, use.spikes=FALSE, subset.row=1:500)
    expect_equal(out3c, trendVar(exprs(X)[setdiff(1:500, which(isSpike(X))),]))
    out3d <- trendVar(X, use.spikes=TRUE, subset.row=1:500)
    expect_equal(out3d, trendVar(exprs(X)[intersect(1:500, which(isSpike(X))),]))
    out3e <- trendVar(X, use.spikes=NA, subset.row=1:500)
    expect_equal(out3e, trendVar(exprs(X)[1:500,]))

    # Checking what happens if all but one feature is a spike-in.
    dummy2 <- rbind(dummy, 0) 
    rownames(dummy2) <- paste0("X", seq_len(nrow(dummy2)))
    X2 <- SingleCellExperiment(list(counts=dummy2))
    isSpike(X2, "Chosen") <- rep(c(TRUE, FALSE), c(ngenes, 1))
    sizeFactors(X2) <- sizeFactors(X2, "Chosen") <- colSums(dummy2)
    X2 <- normalize(X2)
    
    out4 <- trendVar(X2)
    expect_equal(out4$mean, out2$mean)
    expect_equal(out4$var, out2$var)
    expect_equal(out4$trend(test.x), out2$trend(test.x))
    expect_equal(out4$design, out2$design)
})

test_that("trendVar checks size factor centering", {
    isSpike(X, "All") <- !logical(nrow(X))
    subX <- X[,1:10] # Checking that it raises a warning upon subsetting (where the size factors are no longer centered).
    expect_warning(trendVar(subX), "size factors not centred")
    suppressWarnings(rnorm <- normalize(subX))
    expect_warning(trendVar(rnorm), NA)
})

test_that("trendVar works with other trend functions",  {
    mx <- max(out$mean)
    mn <- min(out$mean)

    # Parametric:
    out.semi <- trendVar(d, parametric=TRUE)
    expect_equal(out$mean, out.semi$mean)
    expect_equal(out$var, out.semi$var)
    expect_is(out.semi$trend, "function") 
    expect_true(out.semi$trend(mx) > out.semi$trend(mx+1))
    expect_equal(out.semi$trend(0), 0)
    expect_equal(out.semi$trend(1:10), sapply(1:10, out.semi$trend)) 

    # Spline:
    out.spl <- trendVar(d, method="spline")
    expect_equal(out$mean, out.spl$mean)
    expect_equal(out$var, out.spl$var)
    expect_is(out.spl$trend, "function") 
    expect_equal(out.spl$trend(mx), out.spl$trend(mx+1))
    expect_equal(out.spl$trend(0), 0)
    expect_equal(out.spl$trend(1:10), sapply(1:10, out.spl$trend)) 

    # Results should be the same with/without weighting, as all weights are the same without blocking.
    out.unw <- trendVar(d, weighted=FALSE)
    expect_equal(out.unw$trend(1:100/10), out$trend(1:100/10))
    out.semi.unw <- trendVar(d, weighted=FALSE, parametric=TRUE)
    expect_equal(out.semi.unw$trend(1:100/10), out.semi$trend(1:100/10))
    out.spl.unw <- trendVar(d, weighted=FALSE, method="spline")
    expect_equal(out.spl.unw$trend(1:100/10), out.spl$trend(1:100/10))

    # Checking deprecation warnings and argument specification works (change to default args next time).
    expect_warning(out.sp <- trendVar(d, span=0.2), "deprecated")
    expect_warning(out.sp2 <- trendVar(d, loess.args=list(span=0.2)), NA)
    expect_equal(out.sp$trend(1:50/5), out.sp2$trend(1:50/5))

    expect_warning(out.sp <- trendVar(d, method="spline", df=5), "deprecated")
    expect_warning(out.sp2 <- trendVar(d, method="spline", spline.args=list(df=5)), NA)
    expect_equal(out.sp$trend(1:20/5), out.sp2$trend(1:20/5))
})

set.seed(20002)
test_that("trendVar works with design matrices", {
    design <- matrix(1, ncol(d), 1)
    out2 <- trendVar(d, design=design)
    expect_equal(out2$mean, rowMeans(d))
    expect_equal(out2$var, apply(d, 1, var))
    expect_equal(out2$trend(1:100/10), out$trend(1:100/10))
 
    # A non-trivial design matrix.
    design <- model.matrix(~factor(rep(c(1,2), each=100)))
    out <- trendVar(d, design=design)
    expect_equal(out$mean, rowMeans(d))
    fit <- lm.fit(y=t(d), x=design)
    effects <- fit$effects[-seq_len(ncol(design)),]
    expect_equal(out$var, colMeans(effects^2))
    
    expect_is(out$trend, "function")
    m <- max(out$mean)
    expect_equal(out$trend(m), out$trend(m+1))
    m <- min(out$mean)
    expect_equal(out$trend(0), 0)
    
    expect_equal(out$design, design)
    
    # Trying again with a design matrix with non-trivial pivoting.
    covariate <- 1:200
    design <- model.matrix(~factor(rep(c(1,2), each=100)) + covariate)
    expect_identical(qr(design, LAPACK=TRUE)$pivot, c(3L, 1L, 2L))
    
    out <- trendVar(d, design=design)
    expect_equal(out$mean, rowMeans(d))
    fit <- lm.fit(y=t(d), x=design)
    effects <- fit$effects[-seq_len(ncol(design)),]
    expect_equal(out$var, colMeans(effects^2))

    # Checking that subsetting is properly considered.
    reordering <- sample(nrow(d))
    out3 <- trendVar(d, design=design, subset.row=reordering)
    ref <- trendVar(d[reordering,], design=design)    
    expect_equal(out3$mean, ref$mean)
    expect_equal(out3$var, ref$var)
    expect_equal(out3$trend(0:10/2), ref$trend(0:10/2))
})

set.seed(20003)
test_that("trendVar works with blocking factors", {
    blocking <- sample(3, ncol(d), replace=TRUE)
    fit <- trendVar(d, block=blocking)

    # Checking means are correctly calculated.
    exp.means <- cbind(rowMeans(d[,blocking==1]),
                       rowMeans(d[,blocking==2]),
                       rowMeans(d[,blocking==3]))
    colnames(exp.means) <- 1:3
    expect_equal(exp.means, fit$mean)

    # Checking variances are correctly calculated.
    exp.var <- cbind(apply(d[,blocking==1], 1, var),
                     apply(d[,blocking==2], 1, var),
                     apply(d[,blocking==3], 1, var))
    colnames(exp.var) <- 1:3
    expect_equal(exp.var, fit$var)

    # No real way to check the trend, other than to make sure it works.
    expect_equal(fit$trend(0), 0)
    expect_equal(fit$trend(1:10), sapply(1:10, fit$trend)) # Checking we get consistent results with returned function.
    expect_equal(fit$trend(100:1/20), sapply(100:1/20, fit$trend))  # More checking, reversed order.

    # Checking that blocking with weights actually makes a difference.
    fit.unw <- trendVar(d, block=blocking, weighted=FALSE)
    expect_equal(fit$mean, fit.unw$mean)
    expect_equal(fit$var, fit.unw$var)
    expect_equal(fit$trend(0), fit.unw$trend(0))
    for (i in 1:5) {
        expect_false(isTRUE(all.equal(fit$trend(i), fit.unw$trend(i))))
    }

    # Checking that it correctly ignores blocks without residual d.f.
    blocking[1] <- 0
    refit <- trendVar(d, block=blocking)
    fit.0 <- trendVar(d[,-1], block=blocking[-1])
    expect_equal(refit$mean, cbind(`0`=d[,1], fit.0$mean))
    expect_equal(refit$var, cbind(`0`=NA_real_, fit.0$var))
    expect_equal(refit$trend(0:30/10), fit.0$trend(0:30/10))

    # Checking that subsetting is properly considered.
    reordering <- sample(nrow(d))
    out3 <- trendVar(d, block=blocking, subset.row=reordering)
    ref <- trendVar(d[reordering,], block=blocking)    
    expect_equal(out3$mean, ref$mean)
    expect_equal(out3$var, ref$var)
    expect_equal(out3$trend(0:10/2), ref$trend(0:10/2))
})

# There's a lot of ways it can fail when silly inputs are supplied.

test_that("trendVar throws the correct errors", {
    expect_error(trendVar(d[0,,drop=FALSE], parametric=FALSE), "need at least 2 values for non-parametric curve fitting") # loess fails with empty input vectors.
    expect_error(trendVar(d[0,,drop=FALSE], parametric=TRUE), "need at least 4 values for non-linear curve fitting")

    expect_error(trendVar(d[,0,drop=FALSE]), "no residual d.f. in 'x' for variance estimation") # QR fails straight away.
    expect_error(trendVar(d[,0,drop=FALSE], block=integer(0)), "no residual d.f. in any level of 'block' for variance estimation") # QR fails straight away.
    expect_error(trendVar(d[,0,drop=FALSE], design=cbind(integer(0))), "no residual d.f. in 'design' for variance estimation") # QR fails straight away.
    expect_error(trendVar(d[,1,drop=FALSE], parametric=FALSE), "no residual d.f. in 'x' for variance estimation") # undefined variance with no d.f.

    expect_error(trendVar(d, block=2), "length of 'block'")
    expect_error(trendVar(d, design=cbind(1)), "number of rows in 'design'")
})

####################################################################################################

# Testing the variance decomposition

set.seed(20001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
isSpike(X, "MySpikes") <- rbinom(ngenes, 1, 0.7)==0
sizeFactors(X) <- sizeFactors(X, "MySpikes") <- colSums(dummy)
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
    expect_equal(unname(metadata(out.all)$resid.df), rep(ncells-1L, ngenes))

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
    expect_equal(unname(metadata(out)$resid.df), rep(ncells-1L, sum(!isSpike(X))))

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
    expect_equal(out3$p.value, as.numeric(pnorm(qnorm(pval) %*% resid.df/sum(resid.df))))

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
    expect_identical(testVar(numeric(0), trended, df=df), numeric(0))
    expect_identical(testVar(observed, numeric(0), df=df), numeric(0))
    expect_identical(testVar(observed, trended, df=numeric(0)), rep(NA_real_, length(observed)))
})

####################################################################################################

set.seed(20003)
test_that("combineVar works correctly", {
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

    # Checking averaging of stats.
    N <- c(ncells, ncol(sub.d), ncol(alt.d))
    DF <- c(ncells - 1L, ncol(sub.d) - 3L, ncol(alt.d) - 2L)
    res <- combineVar(dec, dec2, dec3)
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

    # Same results upon subsetting.
    reres <- combineVar(dec[1:10,], dec2[1:10,], dec3[1:10,])
    rescheck <- res[1:10,]
    rescheck$FDR <- p.adjust(rescheck$p.value, method="BH")
    expect_equal(reres, rescheck)

    # Checking fails: 
    dec4 <- dec3
    metadata(dec4) <- list()
    expect_warning(res <- combineVar(dec, dec4), "inputs should come from decomposeVar()", fixed=TRUE)
    expect_error(res <- combineVar(dec, dec2[rev(rownames(dec)),]), "gene identities should be the same") # when you switch up the order.

    redec <- dec
    redec2 <- dec2
    rownames(redec) <- rownames(redec2) <- paste0(rownames(dec), "Y")
    expect_error(res <- combineVar(redec, redec2), "gene names")

    # Just directly returns the input if only one DF is supplied.
    expect_equal(combineVar(dec), dec)
    expect_equal(combineVar(dec2), dec2)
    expect_equal(combineVar(dec3), dec3)
})

####################################################################################################

REFFUN <- function(means, size.factors, tol=1e-6, dispersion=0, pseudo.count=1) {
    # Defining which distribution to use.
    if (dispersion==0) {
        qfun <- function(..., mean) { qpois(..., lambda=mean) }
        dfun <- function(..., mean) { dpois(..., lambda=mean) }
    } else {
        qfun <- function(..., mean) { qnbinom(..., mu=mean, size=1/dispersion) }
        dfun <- function(..., mean) { dnbinom(..., mu=mean, size=1/dispersion) }
    }
    collected.means <- collected.vars <- numeric(length(means))

    for (i in seq_along(means)) {
        m <- means[i]*size.factors
        lower <- qfun(tol, mean=m, lower.tail=TRUE)
        upper <- qfun(tol, mean=m, lower.tail=FALSE)

        # Creating a function to compute the relevant statistics.
        .getValues <- function(j) {
            ranged <- lower[j]:upper[j]
            p <- dfun(ranged, mean=m[j])
            lvals <- log2(ranged/size.factors[j] + pseudo.count)
            return(list(p=p, lvals=lvals))
        }

        # Computing the mean.
        cur.means <- numeric(length(size.factors))
        for (j in seq_along(size.factors)) { 
            out <- .getValues(j)
            cur.means[j] <- sum(out$lvals * out$p) / sum(out$p)
        }
        final.mean <- mean(cur.means)
        collected.means[i] <- final.mean
       
        # Computing the variance. Done separately to avoid
        # storing 'p' and 'lvals' in memory, but as a result
        # we need to compute these values twice.
        cur.vars <- numeric(length(size.factors))
        for (j in seq_along(size.factors)) { 
            out <- .getValues(j)
            cur.vars[j] <- sum((out$lvals - final.mean)^2 * out$p) / sum(out$p)
        }
        collected.vars[i] <- mean(cur.vars)
    }

    return(list(mean=collected.means, var=collected.vars))
}

set.seed(20004)
test_that("makeTechTrend works correctly", {
    # Checking the C++ function against a reference R implementation.
    sf <- runif(100)
    means <- sort(runif(20, 0, 50))
    ref <- REFFUN(means, sf, tol=1e-6, dispersion=0, pseudo.count=1)
    out <- .Call(scran:::cxx_calc_log_count_stats, means, sf, 1e-6, 0, 1)
    expect_equal(ref[[1]], out[[1]])
    expect_equal(ref[[2]], out[[2]])

    ref <- REFFUN(means, sf, tol=1e-6, dispersion=0.1, pseudo.count=1)
    out <- .Call(scran:::cxx_calc_log_count_stats, means, sf, 1e-6, 0.1, 1)
    expect_equal(ref[[1]], out[[1]])
    expect_equal(ref[[2]], out[[2]])
              
    # Testing out all the options.
    out <- makeTechTrend(c(1, 5), pseudo.count=2, size.factors=c(0.5, 1.5))
    log.values <- c(log2(rpois(1e6, lambda=0.5)/0.5 + 2),
                    log2(rpois(1e6, lambda=1.5)/1.5 + 2))
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)
    log.values <- c(log2(rpois(1e6, lambda=0.5*5)/0.5 + 2),
                    log2(rpois(1e6, lambda=1.5*5)/1.5 + 2))
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)

    out <- makeTechTrend(c(1, 5), dispersion=0.1)
    log.values <- log2(rnbinom(1e6, mu=1, size=10) + 1)
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)
    log.values <- log2(rnbinom(1e6, mu=5, size=10) + 1)
    expect_true(abs(out(mean(log.values)) - var(log.values)) < 5e-3)

    # Handles zeroes properly.
    out <- makeTechTrend(0:5)
    expect_equal(out(0), 0)

    # Chucks an error when size factors are not centred.
    expect_error(makeTechTrend(0:5, size.factors=1:5), "centred at unity") 

    # Handles SCE inputs properly.
    X <- SingleCellExperiment(list(counts=matrix(1:5, ncol=2, nrow=5)))
    expect_error(makeTechTrend(x=X), "log.exprs.offset")

    sizeFactors(X) <- c(0.9, 1.1)
    suppressWarnings(X <- normalize(X))
    out <- makeTechTrend(x=X)
    ref <- makeTechTrend(2^seq(0, max(rowMeans(exprs(X))), length.out=100)-1,
                         size.factors=sizeFactors(X))
    expect_equal(out(0:10/2), ref(0:10/2))

    X <- SingleCellExperiment(list(counts=matrix(1:10, ncol=2, nrow=5)))
    suppressWarnings(X <- normalize(X))
    out <- makeTechTrend(x=X)
    libsizes <- colSums(counts(X))
    ref <- makeTechTrend(2^seq(0, max(rowMeans(exprs(X))), length.out=100)-1,
                         size.factors=libsizes/mean(libsizes))
    expect_equal(out(0:10/2), ref(0:10/2))
})

