# This checks the correlateNull function.

# require(scran); require(testthat); source("test-correlate.R")

refnull <- function(niters, ncells) {
    rankings <- as.double(seq_len(ncells))
    shuffled <- .Call(scran:::cxx_auto_shuffle, rankings, as.integer(niters))
    out <- cor(shuffled, rankings, method="spearman")
    sort(out)
}

set.seed(100)
ref <- refnull(1e3, 121)
set.seed(100)
out <- correlateNull(121, iters=1e3)
expect_equal(ref, out)

set.seed(100)
ref <- refnull(1e3, 12)
set.seed(100)
out <- correlateNull(12, iters=1e3)
expect_equal(ref, out)

design <- model.matrix(~factor(rep(c(1,2), each=10)))
df <- nrow(design) - ncol(design)
set.seed(100)
out1 <- correlateNull(df, iters=1e3)
set.seed(100)
out2 <- correlateNull(design=design, iters=1e3)
expect_equal(out1, out2)

# Checking nonsense inputs.

expect_error(correlateNull(ncells=100, iters=0), "number of iterations should be positive")
expect_error(correlateNull(ncells=0), "number of cells should be greater than 2")
expect_error(correlateNull(ncells=100, design=design), "cannot specify both 'ncells' and 'design'")

####################################################################################################
# Checking what happens for the error-tolerant ranking.

whee <- runif(100, -1e-16, 1e-16)
set.seed(100)
r <- scran:::.tolerant_rank(whee)
set.seed(100)
r2 <- rank(integer(100), ties.method="random")
expect_identical(r, r2)

set.seed(200)
extra <- sample(10, 100, replace=TRUE)
set.seed(100)
r <- scran:::.tolerant_rank(whee + extra)
set.seed(100)
r2 <- rank(extra, ties.method="random")
expect_identical(r, r2)

####################################################################################################

checkCorrelations <- function(out, exprs, null.dist) {
    ranked.exprs <- apply(exprs, 1, FUN=rank, ties.method="random")
    colnames(ranked.exprs) <- rownames(exprs)
    assembled.pval <- assembled.rho <- numeric(nrow(out))
    for (p in seq_along(assembled.rho)) { 
        assembled.rho[p] <- cor(ranked.exprs[,out$gene1[p]], ranked.exprs[,out$gene2[p]], method="spearman")        
        assembled.pval[p] <- min(sum(null.dist <= assembled.rho[p] + 1e-8), sum(null.dist >= assembled.rho[p] - 1e-8))
    }
    assembled.pval <- 2*(assembled.pval + 1)/(length(null.dist)+1)
    return(data.frame(rho=assembled.rho, pvalue=assembled.pval, FDR=p.adjust(assembled.pval, method="BH")))
}

####################################################################################################

# This checks the correlatePairs function.

set.seed(10000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

# With a pre-specified distribution.

nulls <- sort(runif(1e6, -1, 1))

set.seed(100)
out <- correlatePairs(X, null.dist=nulls)
set.seed(100)
ref <- checkCorrelations(out, X, null.dist=nulls)

expect_equal(out$rho, ref$rho)
expect_equal(out$p.value, ref$pvalue)
expect_equal(out$FDR, ref$FDR)

# Without a pre-specified distribution.

set.seed(10001)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

set.seed(100)
out <- correlatePairs(X)
set.seed(100)
nulls <- correlateNull(Ncells)
ref <- checkCorrelations(out, X, null.dist=nulls)

expect_equal(out$rho, ref$rho)
expect_equal(out$p.value, ref$pvalue)
expect_equal(out$FDR, ref$FDR)

# With a design matrix.

set.seed(10002)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))
design <- model.matrix(~factor(rep(c(1,2), each=50)))

set.seed(100)
out <- correlatePairs(X, design=design)
set.seed(100)
nulls <- correlateNull(design=design)
fit <- lm.fit(x=design, y=t(X))
exprs <- t(fit$effects[-fit$qr$pivot[seq_len(fit$rank)],])
ref <- checkCorrelations(out, exprs, null.dist=nulls)

expect_equal(out$rho, ref$rho)
expect_equal(out$p.value, ref$pvalue)
expect_equal(out$FDR, ref$FDR)

# Checking that it works with a SCESet object.

set.seed(10003)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))
nulls <- sort(runif(1e6, -1, 1))

set.seed(100)
ref <- correlatePairs(X, null.dist=nulls)
set.seed(100)
X2 <- newSCESet(exprsData=data.frame(X))
out <- correlatePairs(X2, null.dist=nulls)
expect_equal(out, ref)

# With spikes.

isSpike(X2) <- rbinom(Ngenes, 1, 0.6)==0L
set.seed(100)
ref <- correlatePairs(exprs(X2)[!isSpike(X2),,drop=FALSE], null.dist=nulls)
set.seed(100)
out <- correlatePairs(X2, null.dist=nulls)
expect_equal(out, ref)

# Checking nonsense inputs.

expect_error(correlatePairs(X[0,], nulls), "need at least two genes to compute correlations")
expect_error(correlatePairs(X[,0], nulls), "number of cells should be greater than 2")
out <- correlatePairs(X, numeric(0))
expect_equal(out$p.value, rep(1, nrow(out)))
