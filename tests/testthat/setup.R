scramble_vector <- function(x, N, seed=runif(1) * .Machine$integer.max) 
# Iteratively shuffle the vector 'x', using C++'s PRNG.
{
    .Call(scran:::cxx_test_shuffle_vector, x, N, seed)
}

scramble_matrix <- function(x, seed=runif(ncol(x)) * .Machine$integer.max) 
# Shuffle each column of 'x', after setting the seed.
{
    .Call(scran:::cxx_test_shuffle_matrix, x, seed)
}


are_PCs_equal <- function(first, second, tol=1e-8) 
# Check if PCs are equal (other than sign).
{
    expect_identical(dim(first), dim(second))
    relative <- first/second
    expect_true(all(colSums(relative > 0) %in% c(0, nrow(first))))
    expect_true(all(abs(abs(relative)-1) < tol))
}

