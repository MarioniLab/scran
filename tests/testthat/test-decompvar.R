# This tests the decomposeVar() function.
# require(scran); require(testthat); source("test-decompvar.R")

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

test_that("decomposeVar works with a one-way design matrix", {
    g <- factor(rep(c(1,2,3), length.out=ncol(X)))
    fit <- trendVar(X, design=g)
    out <- decomposeVar(X, fit, get.spikes=FALSE)
    expect_equal(out$mean, unname(rowMeans(logcounts(X))[!isSpike(X)]))

    ref.fit <- trendVar(X, design=model.matrix(~g))
    ref.out <- decomposeVar(X, fit, get.spikes=FALSE)
    expect_equal(out, ref.out)

    # Checking what happens when I use design in decomposeVar only.
    fit.none <- trendVar(X)
    out2 <- decomposeVar(X, fit.none, design=g, get.spikes=FALSE)
    expect_equivalent(out2$mean, ref.out$mean)
    expect_equivalent(out2$total, ref.out$total)
    expect_equivalent(out2$tech, fit.none$trend(out$mean))
    expect_equivalent(out2$bio, out2$total-out2$tech)
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


