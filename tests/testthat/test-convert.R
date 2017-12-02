# Tests the convertTo function.
# require(scran); require(testthat); source("test-convert.R")

set.seed(40000)
ncells <- 200
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))

X <- SingleCellExperiment(list(counts=dummy))
is.spike <- rbinom(ngenes, 1, 0.5)==0L
isSpike(X, "MySpike") <- is.spike
sizeFactors(X) <- 2^rnorm(ncells)
expect_warning(X <- normalize(X), "spike-in set 'MySpike'")

rowData(X)$SYMBOL <- paste0("X", seq_len(ngenes))
X$other <- sample(LETTERS, ncells, replace=TRUE)

# Converting to a DGEList.

test_that("Can convert from SingleCellExperiment to DGEList", {
    y <- convertTo(X, type="edgeR")
    expect_identical(y$counts, counts(X)[!is.spike,])
    expect_identical(y$genes, NULL)
    expect_identical(y$offset, NULL) 
    
    y <- convertTo(X, type="edgeR", get.spikes=TRUE)
    expect_identical(y$counts, counts(X))
    expect_identical(y$offset, NULL) # Still NULL, as no spike-in-specific factors are set.
   
    # Checking subsetting behaves as expected. 
    chosen <- c(50:1, 101:200)
    y <- convertTo(X, type="edgeR", subset.row=chosen, get.spikes=TRUE)
    expect_identical(y$counts, counts(X)[chosen,])
    y <- convertTo(X, type="edgeR", subset.row=chosen, get.spikes=FALSE)
    expect_identical(y$counts, counts(X)[setdiff(chosen, which(is.spike)),])
    
    elib <- edgeR::getOffset(y)
    elib <- elib - mean(elib) + mean(log(sizeFactors(X)))
    expect_true(all(abs(exp(elib) - sizeFactors(X)) < 1e-8))
    
    # Checking metadata behaves as expected.
    y <- convertTo(X, type="edgeR", row.fields="SYMBOL", col.fields="other")
    expect_identical(y$samples$other, X$other)
    expect_identical(y$genes$SYMBOL, rowData(X)$SYMBOL[!is.spike])
    
    # Seeing how it behaves with offset matrices.
    X2 <- X 
    sizeFactors(X2, type="MySpike") <- 1
    y <- convertTo(X2, type="edgeR", get.spikes=TRUE)
    log.mls <- mean(log(y$samples$lib.size))
    expect_equal(unname(y$offset[isSpike(X2),,drop=FALSE]), matrix(log.mls, sum(isSpike(X2)), ncol(X2)))
    expect_equal(unname(y$offset[!isSpike(X2),,drop=FALSE]),
                 matrix(log(sizeFactors(X2)) - mean(log(sizeFactors(X2))) + log.mls, 
                        sum(!isSpike(X2)), ncol(X2), byrow=TRUE))
    y <- convertTo(X2, type="edgeR", use.all.sf=FALSE, get.spikes=TRUE)
    expect_identical(y$offset, NULL) 
    y <- convertTo(X2, type="edgeR")
    expect_identical(y$offset, NULL) # NULL because no rows are spike-ins.
   
    # Trying out silly settings. 
    expect_warning(y <- convertTo(X[0,], type="edgeR", row.fields="SYMBOL"))
    expect_identical(y$counts, counts(X)[0,])
    expect_identical(y$genes$SYMBOL, character(0))
    
    y <- convertTo(X[,0], type="edgeR", col.fields="other")
    expect_identical(y$counts, counts(X)[!is.spike,0])
    expect_identical(y$samples$other, character(0))
})

# Converting to a DESeqDataSet.

test_that("Can convert SingleCellExperiment to a DESeqDataSet", {
    library(DESeq2)
    y <- convertTo(X, type="DESeq2")
    expect_equal(counts(y), counts(X)[!is.spike,])
    expect_identical(unname(sizeFactors(y)), sizeFactors(X))
    
    y <- convertTo(X, type="DESeq2", get.spikes=TRUE)
    expect_equal(counts(y), counts(X))
    expect_identical(normalizationFactors(y), NULL) # NULL, as spike-in-specific size factors have not been set.
    
    # Seeing how it behaves with offset matrices.
    X2 <- X 
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
    
    # Checking subsetting behaves as expected.
    chosen <- c(50:1, 101:200)
    y <- convertTo(X, type="DESeq2", subset.row=chosen, get.spikes=TRUE)
    expect_equal(counts(y), counts(X)[chosen,])
    y <- convertTo(X, type="edgeR", subset.row=chosen, get.spikes=FALSE)
    expect_identical(y$counts, counts(X)[setdiff(chosen, which(is.spike)),])

    # Checking metadata is extracted as expected.
    y <- convertTo(X, type="DESeq2", row.fields="SYMBOL", col.fields="other")
    expect_identical(y$other, X$other)
    expect_identical(mcols(y)$SYMBOL, rowData(X)$SYMBOL[!is.spike])

    # # Looks like the DESeqDataSet constructor just fails with no rows/columns.
    # y <- convertTo(X[0,], type="DESeq2", row.fields="SYMBOL")
    # expect_identical(counts(y), counts(X)[0,])
    # expect_identical(mcols(y)$SYMBOL, character(0))
    # 
    # y <- convertTo(X[,0], type="DESeq2", col.fields="other")
    # expect_identical(counts(y), counts(X)[!is.spike,0])
    # expect_identical(y$other, character(0))
})

# Converting to a CellDataSet.

catch_warning <- function(...) {
    expect_warning(..., "gene_short_name")
}

get_exprs <- function(y) {
    assayDataElement(y, "exprs")
}

test_that("Can convert SingleCellExperiment to a CellDataSet", {
    to.comp <- t(t(counts(X))/sizeFactors(X))
    catch_warning(y <- convertTo(X, type="monocle"))
    expect_equal(get_exprs(y), to.comp[!is.spike,])

    # No normalization.
    catch_warning(y <- convertTo(X, type="monocle", normalize=FALSE))
    expect_identical(get_exprs(y), counts(X)[!is.spike,])   
    
    # Assuming no spike-in-specific normalization.
    catch_warning(y <- convertTo(X, type="monocle", get.spikes=TRUE)) 
    expect_equal(get_exprs(y), to.comp)
    
    # Now with spike-in-specific normalization.
    X2 <- X 
    sizeFactors(X2, type="MySpike") <- 1
    catch_warning(y <- convertTo(X2, type="monocle", get.spikes=TRUE))
    expect_equal(get_exprs(y)[!isSpike(X2),], to.comp[!is.spike,])
    expect_equal(get_exprs(y)[isSpike(X2),], counts(X)[is.spike,])
    catch_warning(y <- convertTo(X2, type="monocle", get.spikes=TRUE, use.all.sf=FALSE))
    expect_equal(get_exprs(y), to.comp)

    # Checking that subsetting works as expected.
    chosen <- c(50:1, 101:200)
    catch_warning(y <- convertTo(X, type="monocle", subset.row=chosen, get.spikes=TRUE))
    expect_equal(get_exprs(y), to.comp[chosen,])
    catch_warning(y <- convertTo(X, type="monocle", subset.row=chosen, get.spikes=FALSE))
    expect_equal(get_exprs(y), to.comp[setdiff(chosen, which(is.spike)),])
    
    # Checking metadata is extracted as expected.
    catch_warning(y <- convertTo(X, type="monocle", row.fields="SYMBOL", col.fields="other"))
    expect_identical(y$other, X$other)
    expect_identical(fData(y)$SYMBOL, rowData(X)$SYMBOL[!is.spike])
    
    # # Looks like the CellDataSet constructor just fails with no rows.
    # y <- convertTo(X[0,], type="monocle", row.fields="SYMBOL")
    # expect_identical(get_exprs(y), to.comp[0,])
    # expect_identical(fData(y)$SYMBOL, character(0))
    
    catch_warning(y <- convertTo(X[,0], type="monocle", col.fields="other"))
    expect_identical(get_exprs(y), to.comp[!is.spike,0])
    expect_identical(y$other, character(0))
    
    X2 <- X
    sizeFactors(X2) <- NULL
    expect_error(convertTo(X2, type="monocle"), "size factors not defined for normalization")
})

# Also testing how the methods behave when no spike-ins are specified.

test_that("Conversion succeeds without any spike-ins", { 
    X <- SingleCellExperiment(list(counts=dummy))
    y <- convertTo(X, type="edgeR")
    expect_identical(counts(X), y$counts)
    dds <- convertTo(X, type="DESeq2")
    expect_equal(counts(X), counts(dds))
    sizeFactors(X) <- 1
    catch_warning(y <- convertTo(X, type="monocle"))
    expect_equal(get_exprs(y), counts(X))
})

