#' @importFrom dqrng generateSeedVectors
.setup_pcg_state <- function(per.core) {
    seeds <- streams <- vector("list", length(per.core))
    last <- 0L
    for (i in seq_along(per.core)) { 
        N <- per.core[i]
        seeds[[i]] <- generateSeedVectors(N, nwords=2)
        streams[[i]] <- last + seq_len(N)
        last <- last + N
    }
    list(seeds=seeds, streams=streams)
}
