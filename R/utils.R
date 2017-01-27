.spikeSubset <- function(x, get.spikes) {
    if (!get.spikes) {
        nokeep <- isSpike(x, warning=FALSE)
        if (!is.null(nokeep) && any(nokeep)) {
            return(!nokeep)
        }
    } 
    return(NULL)
}

.subset_to_index <- function(subset, x, byrow=TRUE) {
    if (byrow) {
        dimlen <- nrow(x)
        names <- rownames(x)
    } else {
        dimlen <- ncol(x)
        names <- colnames(x)
    }

    if (is.logical(subset)) { 
        if (length(subset)!=dimlen) {
            stop("subset vector is longer than matrix dimensions") 
        }
        subset <- which(subset)
    } else if (is.character(subset)) {
        subset <- match(subset, names)
        if (any(is.na(subset))) { 
            stop("missing names in subset vector")
        }
    } else if (is.null(subset)) {
        subset <- seq_len(dimlen)
    } else if (is.numeric(subset)) {
        subset <- as.integer(subset)
        if (min(subset) < 1L || max(subset) > dimlen) {
            stop("subset indices out of range")
        }
    } else {
        stop("unrecognized type of subset vector")
    }

    return(subset)
}

.makeVarDefaults <- function(x, fit, design) 
# Makes defaults for the trendVar and decomposeVar functions.
{
    if (is.null(design)) { 
       design <- as.matrix(rep(1, ncol(x))) 
    } else if (length(design)==1L && is.na(design)) { 
        design <- fit$design 
    }
    return(list(design=design))
}

.check_centered_SF <- function(x, assay) 
# Checks if 'exprs' was requested, and if it could have been computed from counts,
# If so, then it checks if the size factors are centered.
{
    if (assay=="exprs" && 
        !is.null(suppressWarnings(get_exprs(x, "counts", warning=TRUE))) && 
        !areSizeFactorsCentred(x)) {
        warning("size factors not centred, run 'normalize()' first")
    }
    return(NULL)
}

.workerAssign <- function(njobs, BPPARAM) 
# Assigns jobs to workers.
{
    ncores <- bpworkers(BPPARAM)
    starting <- as.integer(seq(from=1, to=njobs+1, length.out=ncores+1))
    starting <- unique(starting[seq_len(ncores)])
    ending <- c((starting - 1L)[-1], njobs)
    return(list(start=starting, end=ending))
}

.isOneWay <- function(design) 
# Checks if design matrix is a one-way layout.
{
    if (nrow(design) <= ncol(design)) {
        stop("design matrix has no residual degrees of freedom")
    }
    group <- designAsFactor(design)
    if (nlevels(group) == ncol(design)) {
        # Stripping out groups with only one level.
        groupings <- split(seq_len(nrow(design)), group)
        groupings[lengths(groupings)==1L] <- NULL
        return(groupings)
    } 
    return(NULL)
}

