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
