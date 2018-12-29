scrambler <- function(x, N, reset=TRUE, seed=runif(1) * .Machine$integer.max) 
# A shuffling function using C++'s PRNG.
{
    .Call(scran:::cxx_auto_shuffle, x, N, seed, reset)
}
