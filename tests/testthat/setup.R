all_positive_integers <- function(N) sample(.Machine$integer.max, N, replace=TRUE)

scramble_vector <- function(x, N, seed=all_positive_integers(1L), stream=1) 
# Iteratively shuffle the vector 'x', using C++'s PRNG.
{
    scran:::test_shuffle_vector(x, N, seed, stream)
}

scramble_matrix <- function(x, seed=all_positive_integers(ncol(x)), stream=seq_len(ncol(x))) 
# Shuffle each column of 'x', after setting the seed.
{
    scran:::test_shuffle_matrix(x, seed, stream)
}


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

