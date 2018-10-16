# This tests the various trendVar() options.
# require(scran); require(testthat); source("test-trendvar.R")

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
    isSpike(X, "whee") <- !logical(nrow(X))
    sizeFactors(X, "whee") <- runif(ncol(X))
    expect_warning(trendVar(X), "size factors not centred")
    suppressWarnings(rnorm <- normalize(X))
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

    # Checking that argument specification has some effect.
    out.sp <- trendVar(d, method="loess", loess.args=list(span=0.2))
    out.sp2 <- trendVar(d, method="loess", loess.args=list(span=0.4))
    expect_false(isTRUE(all.equal(out.sp$trend(1:50/5), out.sp2$trend(1:50/5))))

    out.sp <- trendVar(d, method="spline", spline.args=list(df=3))
    out.sp2 <- trendVar(d, method="spline", spline.args=list(df=5))
    expect_false(isTRUE(all.equal(out.sp$trend(1:50/5), out.sp2$trend(1:50/5))))
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
    expect_equal(out$design, design)

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

set.seed(200021)
test_that("trendVar works with one-way design matrices", {
    g <- factor(rep(c(1,2), each=100))
    D <- model.matrix(~g)
    out <- trendVar(d, design=g)
    out2 <- trendVar(d, design=D)

    expect_equal(out$design, g)
    expect_equal(out2$design, D)

    expect_equal(out$vars, out2$vars)
    expect_equal(out$means, out2$means)
    expect_equal(out$trend(0:10/2), out2$trend(0:10/2))
    expect_equal(out$resid.df, out2$resid.df)

    # Checking that subsetting is properly considered.
    reordering <- sample(nrow(d))
    out3 <- trendVar(d, design=g, subset.row=reordering)
    ref <- trendVar(d[reordering,], design=g)
    expect_equal(out3$mean, ref$mean)
    expect_equal(out3$var, ref$var)
    expect_equal(out3$trend(0:10/2), ref$trend(0:10/2))

    # Checking that it ignores levels with no residual d.f.
    g0 <- factor(rep(0:2, c(1, 99, 100)))
    refit <- trendVar(d, design=g0)
    fit.0 <- trendVar(d, design=model.matrix(~g0))
    expect_equal(refit$mean, fit.0$mean)
    expect_equal(refit$var, fit.0$var)
    expect_equal(refit$trend(0:30/10), fit.0$trend(0:30/10))
    expect_equal(refit$resid.df, fit.0$resid.df)
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

test_that("trendVar throws the correct errors", {
    expect_error(trendVar(d[0,,drop=FALSE], parametric=FALSE), "need at least 2 points for non-parametric curve fitting") # loess fails with empty input vectors.
    expect_error(trendVar(d[0,,drop=FALSE], parametric=TRUE), "need at least 4 points for non-linear curve fitting")

    expect_error(trendVar(d[,0,drop=FALSE]), "no residual d.f. in 'x' for variance estimation")
    expect_error(trendVar(d[,0,drop=FALSE], block=integer(0)), "no residual d.f. in any level of 'block' for variance estimation")
    expect_error(trendVar(d[,0,drop=FALSE], design=integer(0)), "no residual d.f. in 'design' for variance estimation")
    expect_error(trendVar(d[,0,drop=FALSE], design=cbind(integer(0))), "no residual d.f. in 'design' for variance estimation")
    expect_error(trendVar(d[,1,drop=FALSE], parametric=FALSE), "no residual d.f. in 'x' for variance estimation")

    expect_error(trendVar(d, block=2), "length of 'block'")
    expect_error(trendVar(d, design=2), "length of one-way 'design'")
    expect_error(trendVar(d, design=cbind(1)), "number of rows in 'design'")
})
