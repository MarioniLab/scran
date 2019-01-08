scramble_vector <- function(x, N, seed=runif(1, 0, 2^32), stream=1) 
# Iteratively shuffle the vector 'x', using C++'s PRNG.
{
    .Call(scran:::cxx_test_shuffle_vector, x, N, seed, stream)
}

scramble_matrix <- function(x, seed=runif(ncol(x), 0, 2^32), stream=seq_len(ncol(x))) 
# Shuffle each column of 'x', after setting the seed.
{
    .Call(scran:::cxx_test_shuffle_matrix, x, seed, stream)
}


are_PCs_equal <- function(first, second, tol=1e-8) 
# Check if PCs are equal (other than sign).
{
    expect_identical(dim(first), dim(second))
    relative <- first/second
    expect_true(all(colSums(relative > 0) %in% c(0, nrow(first))))
    expect_true(all(abs(abs(relative)-1) < tol))
}

