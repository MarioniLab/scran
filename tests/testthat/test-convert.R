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
sizeFactors(X) <- 2^rnorm(ncells)
rowData(X)$SYMBOL <- paste0("X", seq_len(ngenes))
X$other <- sample(LETTERS, ncells, replace=TRUE)

# Converting to a DGEList.

test_that("Can convert from SingleCellExperiment to DGEList", {
    y <- convertTo(X, type="edgeR")
    expect_identical(y$counts, counts(X))
   
    # Checking subsetting behaves as expected. 
    chosen <- c(50:1, 101:200)
    y <- convertTo(X, type="edgeR", subset.row=chosen) 
    expect_identical(y$counts, counts(X)[chosen,])
    
    # Checking metadata behaves as expected.
    y <- convertTo(X, type="edgeR")
    expect_identical(y$samples$other, X$other)
    expect_identical(y$genes$SYMBOL, rowData(X)$SYMBOL)
    
#    # Trying out silly settings. 
#    expect_warning(y <- convertTo(X[0,], type="edgeR"))
#    expect_identical(y$counts, counts(X)[0,])
#    expect_identical(y$genes$SYMBOL, character(0))
#    
#    y <- convertTo(X[,0], type="edgeR")
#    expect_identical(y$counts, counts(X)[,0])
#    expect_identical(y$samples$other, character(0))
})

test_that("Can convert SingleCellExperiment to a DESeqDataSet", {
    library(DESeq2)
    y <- convertTo(X, type="DESeq2")
    expect_equal(counts(y), counts(X))
    expect_identical(unname(sizeFactors(y)), sizeFactors(X))

    # Checking subsetting behaves as expected.
    chosen <- c(50:1, 101:200)
    y <- convertTo(X, type="DESeq2", subset.row=chosen) 
    expect_equivalent(assay(y), counts(X)[chosen,])

    # Checking metadata is extracted as expected.
    y <- convertTo(X, type="DESeq2")
    expect_identical(y$other, X$other)
    expect_identical(mcols(y)$SYMBOL, rowData(X)$SYMBOL)
})

# Converting to a CellDataSet.

catch_warning <- function(...) {
    expect_warning(..., "gene_short_name")
}

get_exprs <- function(y) {
    assayDataElement(y, "exprs")
}

test_that("Can convert SingleCellExperiment to a CellDataSet", {
    skip("monocle is a little broken because of clusterApply")
    catch_warning(y <- convertTo(X, type="monocle"))
    expect_equal(get_exprs(y), counts(X))
    expect_equivalent(sizeFactors(y), sizeFactors(X))

    # Checking that subsetting works as expected.
    chosen <- c(50:1, 101:200)
    catch_warning(y <- convertTo(X, type="monocle", subset.row=chosen))
    expect_equal(get_exprs(y), counts(X)[chosen,])
    expect_equivalent(sizeFactors(y), sizeFactors(X))
    
    # Checking metadata is extracted as expected.
    catch_warning(y <- convertTo(X, type="monocle"))
    expect_identical(y$other, X$other)
    expect_identical(fData(y)$SYMBOL, rowData(X)$SYMBOL)
    
    # # Looks like the CellDataSet constructor just fails with no rows.
    # y <- convertTo(X[0,], type="monocle", row.fields="SYMBOL")
    # expect_identical(get_exprs(y), counts(X)[0,])
    # expect_identical(fData(y)$SYMBOL, character(0))
    
    catch_warning(y <- convertTo(X[,0], type="monocle") )
    expect_identical(get_exprs(y), counts(X)[,0])
    expect_identical(y$other, character(0))
    
    X2 <- X
    sizeFactors(X2) <- NULL
    expect_warning(y2 <- convertTo(X2, type="monocle"))
    expect_equal(counts(X2), get_exprs(y2))
    expect_true(all(is.na(sizeFactors(y2))))
})
