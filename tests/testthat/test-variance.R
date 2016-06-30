# This tests the variance calculation functions in scran.

require(scran); require(testthat);

set.seed(20000)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ngenes*ncells, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

d <- exprs(X)
out <- trendVar(d)
expect_equal(out$mean, rowMeans(d))
expect_equal(out$var, apply(d, 1, var))

expect_is(out$trend, "function") # hard to test it without copying the code, so I'll just check that the bounds are right.
m <- max(out$mean)
expect_equal(out$trend(m), out$trend(m+1))
m <- min(out$mean)
expect_equal(out$trend(m), out$trend(m-1))

expect_equal(out$design, as.matrix(rep(1, ncells)))

# Get the same results directly on a SCESet.

suppressWarnings(expect_error(trendVar(X), "'degree' must be less than number of unique points")) # because there aren't any spike-ins.
isSpike(X) <- TRUE
out2 <- trendVar(X)
expect_equal(out$mean, out2$mean)
expect_equal(out$var, out2$var)
expect_equal(out$trend, out2$trend)
expect_equal(out$design, out2$design)

isSpike(X) <- FALSE
expect_error(trendVar(X), "'degree' must be less than number of unique points")
out3 <- trendVar(X, use.spikes=FALSE)
expect_equal(out3$mean, out2$mean)
expect_equal(out3$var, out2$var)
expect_equal(out3$trend, out2$trend)
expect_equal(out3$design, out2$design)

isSpike(X) <- rbinom(ngenes, 1, 0.5)==0
out3b <- trendVar(X, use.spikes=NA)
expect_equal(out3$mean, out3b$mean)
expect_equal(out3$var, out3b$var)
expect_equal(out3$trend, out3b$trend)
expect_equal(out3$design, out3b$design)

dummy2 <- rbind(dummy, 0)
rownames(dummy2) <- paste0("X", seq_len(nrow(dummy2)))
X2 <- newSCESet(countData=data.frame(dummy2))
sizeFactors(X2) <- colSums(dummy2)
X2 <- normalize(X2)
isSpike(X2) <- rep(c(TRUE, FALSE), c(ngenes, 1))

out4 <- trendVar(X2)
expect_equal(out4$mean, out2$mean)
expect_equal(out4$var, out2$var)
expect_equal(out4$trend, out2$trend)
expect_equal(out4$design, out2$design)

# Trying again but with the loess.

out2 <- trendVar(d, trend="loess")
expect_equal(out$mean, out2$mean)
expect_equal(out$var, out2$var)
expect_equal(out2$design, out$design)

expect_is(out2$trend, "function")
m <- max(out2$mean)
expect_equal(out2$trend(m), out2$trend(m+1))
m <- min(out2$mean)
expect_equal(out2$trend(m), out2$trend(m-1))

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
expect_equal(out$trend(m), out$trend(m-1))

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

# There's a lot of ways it can fail, e.g., trend fitting will not work depending on the number of 'df'.

expect_error(trendVar(d[0,,drop=FALSE]), "'degree' must be less than number of unique points")
expect_error(trendVar(d[2,,drop=FALSE]), "'degree' must be less than number of unique points")
expect_error(trendVar(d[,0,drop=FALSE]), "design matrix is not of full rank")
expect_error(trendVar(d[,1,drop=FALSE]), "design matrix is not of full rank")

####################################################################################################

# Testing the variance decomposition

set.seed(20001)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ngenes*ncells, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- newSCESet(countData=data.frame(dummy))
sizeFactors(X) <- colSums(dummy)
isSpike(X) <- rbinom(ngenes, 1, 0.7)==0
X <- normalize(X)

fit <- trendVar(X)
out <- decomposeVar(X, fit)
out.all <- decomposeVar(X, fit, get.spikes=TRUE)
expect_identical(rownames(out), rownames(X))
expect_identical(rownames(out.all), rownames(X))

ref.mean <- rowMeans(exprs(X))
expect_equivalent(out.all$mean, ref.mean)
ref.mean[isSpike(X)] <- NA_real_
expect_equivalent(out$mean, ref.mean)

ref.var <- apply(exprs(X), 1, var)
expect_true(all(is.na(out.all$total)==is.na(ref.var)))
expect_true(all(abs(unname(out.all$total)-ref.var) < 1e-8 | is.na(ref.var)))
ref.var[isSpike(X)] <- NA_real_
expect_true(all(is.na(out$total)==is.na(ref.var)))
expect_true(all(abs(unname(out$total)-ref.var) < 1e-8 | is.na(ref.var)))

expect_equivalent(out$tech, fit$trend(ref.mean))
expect_equivalent(out$bio, out$total-out$tech)
expect_true(all(abs(out.all$tech - fit$trend(out.all$mean)) < 1e-8 | is.na(out.all$tech)))
expect_equivalent(out.all$bio, out.all$total-out.all$tech)

shuffled <- c(500:1, 501:1000)
out.ref <- decomposeVar(X[shuffled,], fit)
out2 <- decomposeVar(X, fit, subset.row=shuffled)
expect_identical(out.ref, out2) # Checking that subset.row works. 

# Testing with a modified design matrix.

out2 <- decomposeVar(X, fit, design=NULL) # defaults to all-ones.
expect_equal(out, out2)

design <- model.matrix(~factor(rep(c(1,2), each=100)))
out3 <- decomposeVar(X, fit, design=design)
expect_equal(out$mean, out3$mean)

refit <- lm.fit(y=t(exprs(X)), x=design)
effects <- refit$effects[-seq_len(ncol(design)),]
test.var <- colMeans(effects^2)
test.var[isSpike(X)] <- NA
expect_equivalent(out3$total, test.var)

expect_equivalent(out3$tech, fit$trend(ref.mean))
expect_equivalent(out3$bio, out3$total-out3$tech)

####################################################################################################

# Testing the testVar() function.

set.seed(20002)
true.p <- runif(100)
trended <- runif(100, 1, 2)
df <- 20
observed <- trended * qchisq(true.p, df=df, lower.tail=FALSE)/df
pvals <- testVar(observed, trended, df=df, min=0)
expect_equal(pvals, true.p)

design <- model.matrix(~factor(rep(c(1,2), each=11)))
pvals <- testVar(observed, trended, design=design, min=0)
expect_equal(pvals, true.p)

# Checking silly inputs

expect_identical(testVar(numeric(0), trended, df=df), numeric(0))
expect_identical(testVar(observed, numeric(0), df=df), rep(NA_real_, length(observed)))
expect_identical(testVar(observed, trended, df=numeric(0)), rep(NA_real_, length(observed)))

