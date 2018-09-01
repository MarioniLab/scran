.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3) 
# Centralized function to compute tricube averages.
# Bandwidth is set at 'ndist' times the median distance, if not specified.
{
    if (is.null(bandwidth)) {
        middle <- ceiling(ncol(indices)/2L)
        mid.dist <- distances[,middle]
        bandwidth <- mid.dist * ndist
    }
    bandwidth <- pmax(1e-8, bandwidth)

    rel.dist <- distances/bandwidth
    rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
    tricube <- (1 - rel.dist^3)^3
    weight <- tricube/rowSums(tricube)

    output <- 0
    for (kdx in seq_len(ncol(indices))) {
        output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
    }

    if (is.null(dim(output))) {
        matrix(0, nrow(vals), ncol(vals))
    } else {
        output
    }
}
