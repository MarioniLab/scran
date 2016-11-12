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

suppressWarnings(expect_error(trendVar(X), "invalid 'x'")) # because there aren't any spike-ins for loess.

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
expect_error(trendVar(X), "invalid 'x'")
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

dummy2 <- rbind(dummy, 0)
rownames(dummy2) <- paste0("X", seq_len(nrow(dummy2)))
X2 <- newSCESet(countData=data.frame(dummy2))
X2 <- calculateQCMetrics(X2, list(Chosen=rep(c(TRUE, FALSE), c(ngenes, 1))))
setSpike(X2) <- "Chosen"
sizeFactors(X2) <- colSums(dummy2)
X2 <- normalize(X2)

out4 <- trendVar(X2)
expect_equal(out4$mean, out2$mean)
expect_equal(out4$var, out2$var)
expect_equal(out4$trend, out2$trend)
expect_equal(out4$design, out2$design)

# Trying again but with the semiloessnomial.

suppressWarnings(out2 <- trendVar(d, trend="semiloess"))
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

# There's a lot of ways it can fail when silly inputs are supplied.

expect_error(trendVar(d[0,,drop=FALSE]), "invalid 'x'") # loess fails with empty input vectors.
expect_error(trendVar(d[0,,drop=FALSE], trend="semiloess"), "need at least 4 values for non-linear curve fitting")
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
X <- calculateQCMetrics(X, list(MySpikes=rbinom(ngenes, 1, 0.7)==0))
setSpike(X) <- "MySpikes"
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

fit <- trendVar(X)
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

expect_equivalent(out3$total, test.var)
expect_equivalent(out3$tech, fit$trend(ref.mean))
expect_equivalent(out3$bio, out3$total-out3$tech)

ref.p <- testVar(out3$total, out3$tech, df=nrow(design) - ncol(design))
ref.p[isSpike(X)] <- NA_real_
expect_equivalent(ref.p, out3$p.value)
expect_equivalent(p.adjust(ref.p, method="BH"), out3$FDR)

####################################################################################################

# Testing the testVar() function.

set.seed(20002)
true.p <- runif(100)
trended <- runif(100, 1, 2)
df <- 20
observed <- trended * qchisq(true.p, df=df, lower.tail=FALSE)/df
pvals <- testVar(observed, trended, df=df)
expect_equal(pvals, true.p)

design <- model.matrix(~factor(rep(c(1,2), each=11)))
pvals <- testVar(observed, trended, design=design)
expect_equal(pvals, true.p)

# Checking silly inputs

expect_identical(testVar(numeric(0), trended, df=df), numeric(0))
expect_identical(testVar(observed, numeric(0), df=df), rep(NA_real_, length(observed)))
expect_identical(testVar(observed, trended, df=numeric(0)), rep(NA_real_, length(observed)))

