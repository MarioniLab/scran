.compute_tricube_average <- function(vals, indices, distances, bandwidth) 
# Centralized function to compute tricube averages.
{
    rel.dist <- pmin(1, distances/bandwidth)
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
