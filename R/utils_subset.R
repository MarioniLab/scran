.subset_to_index <- function(subset, x, byrow=TRUE) 
# Converts arbitrary subsetting vectors to an integer index vector.
{
    if (byrow) {
        dummy <- seq_len(nrow(x))
        names(dummy) <- rownames(x)
    } else {
        dummy <- seq_len(ncol(x))
        names(dummy) <- colnames(x) 
    }

    if (!is.null(subset)) { 
        dummy <- dummy[subset]
    }
    out <- unname(dummy)

    if (any(is.na(out))) {
        stop("'subset' indices out of range of 'x'")
    }
    out
}
