#' @importFrom BiocParallel bpnworkers
.worker_assign <- function(njobs, BPPARAM)
# Assigns jobs to workers, where the each element of the output list is 
# the index vector of the jobs. These are guaranteed to be consecutive
# so the bplapply's output can just be directly combined.
{
    ncores <- bpnworkers(BPPARAM)
    starting <- as.integer(seq(from=1, to=njobs+1, length.out=ncores+1))
    jobsize <- diff(starting)
    starting <- starting[-length(starting)] - 1L
    return(mapply("+", starting, lapply(jobsize, seq_len), SIMPLIFY=FALSE))
}

.split_vector_by_workers <- function(vec, assignments) 
# Convenience function to split a vector by the assigned core,
# to get something that can be easily used in 'bplapply'.
{
    by.core <- vector("list", length(assignments))
    for (core in seq_along(assignments)) {
        by.core[[core]] <- vec[assignments[[core]]]
    }
    names(by.core) <- names(assignments)
    return(by.core)
}

.split_matrix_by_workers <- function(mat, assignments, byrow=TRUE) 
# Convenience function to split a mat by the assigned core,
# to get something that can be easily used in 'bplapply'.
{
    by.core <- vector("list", length(assignments))
    for (core in seq_along(assignments)) {
        chosen <- assignments[[core]]
        if (byrow) {
            current <- mat[chosen,,drop=FALSE]
        } else {
            current <- mat[,chosen,drop=FALSE]
        }
        by.core[[core]] <- current
    }
    names(by.core) <- names(assignments)
    return(by.core)
}
