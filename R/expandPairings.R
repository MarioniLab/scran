.expand_pairings <- function(pairings, universe) {
    .SUBSET <- function(request, clean=TRUE) {
        if (is.null(request)) {
            out <- seq_along(universe)
        } else { 
            out <- match(request, universe)
        }
        if (clean) {
            out <- unique(out[!is.na(out)])
        }
        out
    }

    .expand_pairings_core(pairings, .SUBSET)
}

.expand_pairings_core <- function(pairings, .SUBSET) {
    .clean_expand <- function(x, y, keep.perm) {
        all.pairs <- expand.grid(x, y)
        keep <- all.pairs[,1] != all.pairs[,2]
        all.pairs[keep,]
    }

    if (is.matrix(pairings)) {
        # If matrix, we're using pre-specified pairs.
        if ((!is.numeric(pairings) && !is.character(pairings)) || ncol(pairings)!=2L) {
            stop("'pairings' should be a numeric/character matrix with 2 columns")
        }
        s1 <- .SUBSET(pairings[,1], clean=FALSE)
        s2 <- .SUBSET(pairings[,2], clean=FALSE)

        # Discarding pairs with missing elements silently.
        keep <- !is.na(s1) & !is.na(s2)
        s1 <- s1[keep]
        s2 <- s2[keep]
        mode <- "predefined pairs"

    } else if (is.list(pairings)) {
        # If list, we're correlating between one gene selected from each of two pools.
        if (length(pairings)!=2L) {
            stop("'pairings' as a list should have length 2")
        }
        converted <- lapply(pairings, FUN=.SUBSET)
        all.pairs <- .clean_expand(converted[[1]], converted[[2]])
        s1 <- all.pairs[,1]
        s2 <- all.pairs[,2]
        mode <- "double pool"

    } else {
        available <- .SUBSET(pairings)
        all.pairs <- .clean_expand(available, available)
        s1 <- all.pairs[,1]
        s2 <- all.pairs[,2]
        mode <- "single pool"
    }

    list(id1=s1, id2=s2, mode=mode)
}
