all_positive_integers <- function(N) sample(.Machine$integer.max, N, replace=TRUE)

are_PCs_equal <- function(first, second, tol=1e-8) 
# Check if PCs are equal (other than sign).
{
    expect_identical(dim(first), dim(second))
    relative <- first/second
    expect_true(all(colSums(relative > 0) %in% c(0, nrow(first))))
    expect_true(all(abs(abs(relative)-1) < tol))
}

# Because SnowParam() is too slow, yet MulticoreParam() fails on Windows.
# See discussion at https://github.com/Bioconductor/BiocParallel/issues/98.
safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}

# Using ExactParam to avoid the trouble of setting the seed for all SVD-related tests.
options(BiocSingularParam.default=BiocSingular::ExactParam())

# Adding a test to flush out any uncontrolled parallelization.
library(BiocParallel)
failgen <- setRefClass("FailParam", 
    contains="BiocParallelParam",     
    fields=list(),
    methods=list())

FAIL <- failgen()
register(FAIL) 

library(DelayedArray)
setAutoBPPARAM(FAIL)
