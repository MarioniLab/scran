#' @export
cosineNorm <- function(X, mode=c("matrix", "all", "l2norm"))
# Computes the cosine norm, with some protection from zero-length norms.
#
# written by Aaron Lun
# 5 July 2018
{
    mode <- match.arg(mode)
    out <- .Call(cxx_cosine_norm, X, mode!="l2norm")
    names(out) <- c("matrix", "l2norm")
    switch(mode, all=out, matrix=out$matrix, l2norm=out$l2norm)
}
