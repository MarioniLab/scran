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
X <- calculateQCMetrics(X, list(MySpike=is.spike))
setSpike(X) <- "MySpike"
sizeFactors(X) <- 2^rnorm(ncells)
expect_warning(X <- normalize(X), "spike-in transcripts in 'MySpike'")

fData(X)$SYMBOL <- paste0("X", seq_len(ngenes))
X$other <- sample(LETTERS, ncells, replace=TRUE)

# Converting to a DGEList.

y <- convertTo(X, type="edgeR")
expect_identical(y$counts, counts(X)[!is.spike,])
expect_identical(y$genes, NULL)
expect_identical(y$offset, NULL) 

y <- convertTo(X, type="edgeR", get.spikes=TRUE)
expect_identical(y$counts, counts(X))
expect_identical(y$offset, NULL) # Still NULL, as no spike-in-specific factors are set.

chosen <- c(50:1, 101:200)
y <- convertTo(X, type="edgeR", subset.row=chosen)
expect_identical(y$counts, counts(X)[chosen,])

elib <- edgeR::getOffset(y)
elib <- elib - mean(elib) + mean(log(sizeFactors(X)))
expect_true(all(abs(exp(elib) - sizeFactors(X)) < 1e-8))

y <- convertTo(X, type="edgeR", fData.col="SYMBOL", pData.col="other")
expect_identical(y$samples$other, X$other)
expect_identical(y$genes$SYMBOL, fData(X)$SYMBOL[!is.spike])

X2 <- X # Seeing how it behaves with offset matrices.
sizeFactors(X2, type="MySpike") <- 1
y <- convertTo(X2, type="edgeR", get.spikes=TRUE)
expect_equal(unname(y$offset[isSpike(X2),,drop=FALSE]), matrix(0, sum(isSpike(X2)), ncol(X2)))
expect_equal(unname(y$offset[!isSpike(X2),,drop=FALSE]),
             matrix(log(sizeFactors(X2)) - mean(log(sizeFactors(X2))), 
                    sum(!isSpike(X2)), ncol(X2), byrow=TRUE))
y <- convertTo(X2, type="edgeR", use.all.sf=FALSE, get.spikes=TRUE)
expect_identical(y$offset, NULL) 
y <- convertTo(X2, type="edgeR")
expect_identical(y$offset, NULL) # NULL because no rows are spike-ins.

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
expect_identical(normalizationFactors(y), NULL) # NULL, as spike-in-specific size factors have not been set.

X2 <- X # Seeing how it behaves with offset matrices.
sizeFactors(X2, type="MySpike") <- 1
y <- convertTo(X2, type="DESeq2", get.spikes=TRUE)
expect_equal(unname(normalizationFactors(y)[isSpike(X2),,drop=FALSE]), matrix(1, sum(isSpike(X2)), ncol(X2)))
expect_equal(unname(normalizationFactors(y)[!isSpike(X2),,drop=FALSE]),
             matrix(sizeFactors(X2)/exp(mean(log(sizeFactors(X2)))), 
                    sum(!isSpike(X2)), ncol(X2), byrow=TRUE))
y <- convertTo(X2, type="DESeq2", use.all.sf=FALSE, get.spikes=TRUE)
expect_identical(normalizationFactors(y), NULL)
y <- convertTo(X2, type="DESeq2")
expect_identical(normalizationFactors(y), NULL)

chosen <- c(50:1, 101:200)
y <- convertTo(X, type="DESeq2", subset.row=chosen)
expect_equal(counts(y), counts(X)[chosen,])

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

y <- convertTo(X, type="monocle", get.spikes=TRUE) # Assuming no spike-in-specific normalization.
expect_equal(exprs(y), to.comp)

X2 <- X # Now, with spike-in-specific normalization.
sizeFactors(X2, type="MySpike") <- 1
y <- convertTo(X2, type="monocle", get.spikes=TRUE)
expect_equal(exprs(y)[!isSpike(X2),], to.comp[!is.spike,])
expect_equal(exprs(y)[isSpike(X2),], counts(X)[is.spike,])
y <- convertTo(X2, type="monocle", get.spikes=TRUE, use.all.sf=FALSE)
expect_equal(exprs(y), to.comp)

chosen <- c(50:1, 101:200)
y <- convertTo(X, type="monocle", subset.row=chosen)
expect_equal(exprs(y), to.comp[chosen,])

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

# Also testing how the methods behave when no spike-ins are specified.

X <- newSCESet(countData=data.frame(dummy))
y <- convertTo(X, type="edgeR")
expect_identical(counts(X), y$counts)
dds <- convertTo(X, type="DESeq2")
expect_equal(counts(X), counts(dds))
sizeFactors(X) <- 1
y <- convertTo(X, type="monocle")
expect_equal(exprs(y), counts(X))

