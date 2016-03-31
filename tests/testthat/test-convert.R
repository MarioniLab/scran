# Tests the convertTo function.

# require(scran); require(testthat)

set.seed(40000)
ncells <- 200
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- newSCESet(countData=data.frame(dummy))
is.spike <- rbinom(ngenes, 1, 0.5)==0L
isSpike(X) <- is.spike
sizeFactors(X) <- 2^rnorm(ncells)
X <- normalize(X)

fData(X)$SYMBOL <- paste0("X", seq_len(ngenes))
X$other <- sample(LETTERS, ncells, replace=TRUE)

# Converting to a DGEList.

y <- convertTo(X, type="edgeR")
expect_identical(y$counts, counts(X)[!is.spike,])
expect_identical(y$genes, NULL)

y <- convertTo(X, type="edgeR", get.spikes=TRUE)
expect_identical(y$counts, counts(X))

elib <- edgeR::getOffset(y)
elib <- elib - mean(elib) + mean(log(sizeFactors(X)))
expect_true(all(abs(exp(elib) - sizeFactors(X)) < 1e-8))

y <- convertTo(X, type="edgeR", fData.col="SYMBOL", pData.col="other")
expect_identical(y$samples$other, X$other)
expect_identical(y$genes$SYMBOL, fData(X)$SYMBOL[!is.spike])

expect_warning(y <- convertTo(X[0,], type="edgeR", fData.col="SYMBOL"))
expect_identical(y$counts, counts(X)[0,])
expect_identical(y$genes$SYMBOL, character(0))

y <- convertTo(X[,0], type="edgeR", pData.col="other")
expect_identical(y$counts, counts(X)[!is.spike,0])
expect_identical(y$samples$other, character(0))

# Converting to a DESeqDataSet.

library(DESeq2)
y <- convertTo(X, type="DESeq2")
expect_equal(counts(y), counts(X)[!is.spike,])
expect_identical(sizeFactors(y), sizeFactors(X))

y <- convertTo(X, type="DESeq2", get.spikes=TRUE)
expect_equal(counts(y), counts(X))

y <- convertTo(X, type="DESeq2", fData.col="SYMBOL", pData.col="other")
expect_identical(y$other, X$other)
expect_identical(mcols(y)$SYMBOL, fData(X)$SYMBOL[!is.spike])

# # Looks like the DESeqDataSet constructor just fails with no rows/columns.
# y <- convertTo(X[0,], type="DESeq2", fData.col="SYMBOL")
# expect_identical(counts(y), counts(X)[0,])
# expect_identical(mcols(y)$SYMBOL, character(0))
# 
# y <- convertTo(X[,0], type="DESeq2", pData.col="other")
# expect_identical(counts(y), counts(X)[!is.spike,0])
# expect_identical(y$other, character(0))

# Converting to a CellDataSet.

to.comp <- t(t(counts(X))/sizeFactors(X)) 
y <- convertTo(X, type="monocle")
expect_equal(exprs(y), to.comp[!is.spike,])

y <- convertTo(X, type="monocle", get.spikes=TRUE)
expect_equal(exprs(y), to.comp)

y <- convertTo(X, type="monocle", normalize=FALSE)
expect_identical(exprs(y), counts(X)[!is.spike,])

y <- convertTo(X, type="monocle", fData.col="SYMBOL", pData.col="other")
expect_identical(y$other, X$other)
expect_identical(fData(y)$SYMBOL, fData(X)$SYMBOL[!is.spike])

# # Looks like the CellDataSet constructor just fails with no rows.
# y <- convertTo(X[0,], type="monocle", fData.col="SYMBOL")
# expect_identical(exprs(y), to.comp[0,])
# expect_identical(fData(y)$SYMBOL, character(0))

y <- convertTo(X[,0], type="monocle", pData.col="other")
expect_identical(exprs(y), to.comp[!is.spike,0])
expect_identical(y$other, character(0))

X2 <- X
sizeFactors(X2) <- NULL
expect_error(convertTo(X2, type="monocle"), "size factors not defined for normalization")

