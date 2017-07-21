# This tests the variance calculation functions in scran.
# require(scran); require(testthat); source("test-variance.R")

set.seed(20000)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)
d <- exprs(X)
out <- trendVar(d)

test_that("trendVar works on a basic scenario", {
    expect_equal(out$mean, rowMeans(d))
    expect_equal(out$var, apply(d, 1, var))
    
    # hard to test it without copying the code, so I'll just check that the bounds are right.
    expect_is(out$trend, "function") 
    mx <- max(out$mean)
    mn <- min(out$mean)
    expect_equal(out$trend(mx), out$trend(mx+1))
    expect_equal(out$trend(0), 0)
    
    expect_equal(out$design, as.matrix(rep(1, ncells)))
})

test_that("trendVar is robust to zero-variance genes", {
    # Checking that genes with no variance don't break the results, but still get reported in the output.
    dz <- rbind(1, d, 0)
    outz <- trendVar(dz)
    expect_equal(outz$mean, c(1, out$mean, 0))
    expect_equal(outz$var, c(0, out$var, 0))
    expect_equal(out$trend(outz$mean), outz$trend(outz$mean))
})

test_that("trendVar behaves correctly with subsetting", {
    shuffled <- c(500:110)
    out.ref <- trendVar(d[shuffled,])
    out2 <- trendVar(d, subset.row=shuffled)
    expect_equal(out.ref, out2) # Checking that subset.row works. 
})

test_that("trendVar works correctly on SCESet objects", {
    # Get the same results directly on a SCESet, with various spike-in specifications.
    suppressWarnings(expect_error(trendVar(X, parametric=TRUE), "need at least 4 values for non-linear curve fitting")) 
    suppressWarnings(expect_error(trendVar(X, parametric=FALSE), "need at least 2 values for non-parametric curve fitting"))

    cntrl_data <- list(All=!logical(ngenes), Some=rbinom(ngenes, 1, 0.5)==0, None=logical(ngenes))
    X <- calculateQCMetrics(X, cntrl_data)
    setSpike(X) <- "All"
    expect_identical(isSpike(X), cntrl_data$All) # Just checking here...
    
    out2 <- trendVar(X)
    expect_equal(out$mean, out2$mean)
    expect_equal(out$var, out2$var)
    expect_equal(out$trend, out2$trend)
    expect_equal(out$design, out2$design)
    
    setSpike(X) <- "None"
    expect_identical(isSpike(X), cntrl_data$None)
    expect_error(trendVar(X), "need at least 2 values for non-parametric curve fitting")
    out3 <- trendVar(X, use.spikes=FALSE)
    expect_equal(out3$mean, out2$mean)
    expect_equal(out3$var, out2$var)
    expect_equal(out3$trend, out2$trend)
    expect_equal(out3$design, out2$design)
    
    setSpike(X) <- "Some"
    expect_identical(isSpike(X), cntrl_data$Some)
    out3a <- trendVar(X)
    expect_equal(out3$mean[cntrl_data$Some], out3a$mean)
    expect_equal(out3$var[cntrl_data$Some], out3a$var)
    out3b <- trendVar(X, use.spikes=NA)
    expect_equal(out3$mean, out3b$mean)
    expect_equal(out3$var, out3b$var)
    expect_equal(out3$trend, out3b$trend)
    expect_equal(out3$design, out3b$design)
    
    dummy2 <- rbind(dummy, 0) # Checking what happens if all but one feature is a spike-in.
    rownames(dummy2) <- paste0("X", seq_len(nrow(dummy2)))
    X2 <- newSCESet(countData=data.frame(dummy2))
    X2 <- calculateQCMetrics(X2, list(Chosen=rep(c(TRUE, FALSE), c(ngenes, 1))))
    setSpike(X2) <- "Chosen"
    sizeFactors(X2) <- sizeFactors(X2, "Chosen") <- colSums(dummy2)
    X2 <- normalize(X2)
    
    out4 <- trendVar(X2)
    expect_equal(out4$mean, out2$mean)
    expect_equal(out4$var, out2$var)
    expect_equal(out4$trend, out2$trend)
    expect_equal(out4$design, out2$design)
})

test_that("trendVar checks size factor centering", {
    X <- calculateQCMetrics(X, list(All=!logical(nrow(X))))
    setSpike(X) <- "All"
    subX <- X[,1:10] # Checking that it raises a warning upon subsetting (where the size factors are no longer centered).
    expect_warning(trendVar(subX), "size factors not centred")
    suppressWarnings(rnorm <- normalize(subX))
    expect_warning(trendVar(rnorm), NA)
})

test_that("trendVar works with other trend functions",  {
    mx <- max(out$mean)
    mn <- min(out$mean)

    out.semi <- trendVar(d, parametric=TRUE)
    expect_equal(out$mean, out.semi$mean)
    expect_equal(out$var, out.semi$var)
    expect_is(out.semi$trend, "function") 
    expect_true(out.semi$trend(mx) > out.semi$trend(mx+1))
    expect_equal(out.semi$trend(0), 0)

    out.spl <- trendVar(d, method="spline")
    expect_equal(out$mean, out.spl$mean)
    expect_equal(out$var, out.spl$var)
    expect_is(out.spl$trend, "function") 
    expect_equal(out.spl$trend(mx), out.spl$trend(mx+1))
    expect_equal(out.spl$trend(0), 0)
})

# Trying again with a design matrix.

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

# There's a lot of ways it can fail when silly inputs are supplied.

test_that("trendVar throws the correct errors", {
    expect_error(trendVar(d[0,,drop=FALSE], parametric=FALSE), "need at least 2 values for non-parametric curve fitting") # loess fails with empty input vectors.
    expect_error(trendVar(d[0,,drop=FALSE], parametric=TRUE), "need at least 4 values for non-linear curve fitting")
    expect_error(trendVar(d[,0,drop=FALSE]), "BLAS/LAPACK routine 'DGEQP3' gave error code -4") # QR fails straight away.
    expect_error(trendVar(d[,1,drop=FALSE], parametric=FALSE), "invalid 'x'") # undefined variance with no d.f.
})

####################################################################################################

# Testing the variance decomposition

set.seed(20001)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- newSCESet(countData=data.frame(dummy))
X <- calculateQCMetrics(X, list(MySpikes=rbinom(ngenes, 1, 0.7)==0))
setSpike(X) <- "MySpikes"
sizeFactors(X) <- sizeFactors(X, "MySpikes") <- colSums(dummy)
X <- normalize(X)

fit <- trendVar(X)

test_that("Variance decomposition is working correctly", {
    out <- decomposeVar(X, fit)
    out.all <- decomposeVar(X, fit, get.spikes=TRUE)
    expect_identical(rownames(out), rownames(X))
    expect_identical(rownames(out.all), rownames(X))
    basic <- c("mean", "total", "bio", "tech")
    expect_equivalent(out.all[,basic], out[,basic]) 
    
    ref.mean <- rowMeans(exprs(X))
    expect_equivalent(out.all$mean, ref.mean)
    ref.var <- apply(exprs(X), 1, var)
    expect_equivalent(ref.var, out.all$total)
    expect_equivalent(out$tech, fit$trend(ref.mean))
    expect_equivalent(out$bio, out$total-out$tech)
    
    ref.p <- testVar(out$total, out$tech, df=ncells-1)
    expect_equivalent(ref.p, out.all$p.value)
    expect_equivalent(p.adjust(ref.p, method="BH"), out.all$FDR)
    ref.p[isSpike(X)] <- NA_real_
    expect_equivalent(ref.p, out$p.value)
    expect_equivalent(p.adjust(ref.p, method="BH"), out$FDR)
})
   
test_that("decomposeVar behaves correctly with subsetting", {
    shuffled <- c(500:1, 501:1000)
    out.ref <- decomposeVar(X[shuffled,], fit)
    out2 <- decomposeVar(X, fit, subset.row=shuffled)
    expect_identical(out.ref, out2) # Checking that subset.row works. 
})

test_that("decomposeVar checks size factor centering", {
    subX <- X[,1:10] # Checking that it raises a warning upon subsetting (where the size factors are no longer centered).
    expect_warning(decomposeVar(subX, fit, design=NULL), "size factors not centred")
    expect_warning(decomposeVar(normalize(subX), fit, design=NULL), NA)
})

test_that("decomposeVar works with all genes", {
    # Using all genes for trend fitting.
    all.fit <- trendVar(X, use.spikes=NA)
    all.dec <- decomposeVar(X, all.fit, get.spikes=TRUE)
    expect_equal(all.fit$mean, setNames(all.dec$mean, rownames(all.dec)))
    expect_equal(all.fit$var, setNames(all.dec$total, rownames(all.dec)))
    expect_equal(all.fit$trend(all.fit$mean), setNames(all.dec$tech, rownames(all.dec)))

    all.dec2 <- decomposeVar(fit=all.fit)
    expect_equal(all.dec, all.dec2)
})

test_that("decomposeVar works with design matrices", {
    # Testing with a modified design matrix.
    fit <- trendVar(X)
    out <- decomposeVar(X, fit)
    out2 <- decomposeVar(X, fit, design=NULL) # defaults to all-ones.
    expect_equal(out, out2)

    design <- model.matrix(~factor(rep(c(1,2), each=100)))
    out3 <- decomposeVar(X, fit, design=design)
    expect_equal(out$mean, out3$mean)

    refit <- lm.fit(y=t(exprs(X)), x=design)
    effects <- refit$effects[-seq_len(ncol(design)),]
    test.var <- colMeans(effects^2)
    
    expect_equivalent(out3$total, test.var)
    expect_equivalent(out3$tech, fit$trend(out$mean))
    expect_equivalent(out3$bio, out3$total-out3$tech)
    
    ref.p <- testVar(out3$total, out3$tech, df=nrow(design) - ncol(design))
    ref.p[isSpike(X)] <- NA_real_
    expect_equivalent(ref.p, out3$p.value)
    expect_equivalent(p.adjust(ref.p, method="BH"), out3$FDR)
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
    df1 <- nrow(fit$design)-ncol(fit$design)
    ffit <- limma::fitFDistRobustly(rat[fit$var > 0 & fit$mean > 0.1], df=df1) # filtering out zero-variance/low-abundance genes.
    expect_equal(pvals, pf(rat/ffit$scale, df1=df1, df2=ffit$df2, lower.tail=FALSE))
    
    expect_error(testVar(fit$var, fit$trend(fit$mean), df=ncells-1, test='f'),
                 "second df from trendVar() must be specified for test='f'", fixed=TRUE)
})

test_that("testVar works with silly inputs", {
    expect_identical(testVar(numeric(0), trended, df=df), numeric(0))
    expect_identical(testVar(observed, numeric(0), df=df), rep(NA_real_, length(observed)))
    expect_identical(testVar(observed, trended, df=numeric(0)), rep(NA_real_, length(observed)))
})

