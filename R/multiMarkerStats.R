#' Combine multiple sets of marker statistics
#'
#' Combine multiple sets of marker statistics, typically from different tests,
#' into a single \linkS4class{DataFrame} for convenient inspection.
#'
#' @param ... Two or more lists or \linkS4class{List}s produced by \code{\link{findMarkers}} or \code{\link{combineMarkers}}.
#' Each list should contain \linkS4class{DataFrame}s of results, one for each group/cluster of cells.
#'
#' The names of each List should be the same; the universe of genes in each DataFrame should be the same;
#' and the same number of columns in each DataFrame should be named.
#' All elements in \code{...} are also expected to be named.
#' @param repeated Character vector of columns that are present in one or more DataFrames but should only be reported once.
#' Typically used to avoid reporting redundant copies of annotation-related columns.
#' @param sorted Logical scalar indicating whether each output DataFrame should be sorted by some relevant statistic.
#'
#' @return
#' A named List of DataFrames with one DataFrame per group/cluster.
#' Each DataFrame contains statistics from the corresponding entry of each List in \code{...},
#' prefixed with the name of the List.
#' In addition, several combined statistics are reported:
#' \itemize{
#' \item \code{Top}, the largest rank of each gene across all DataFrames for that group.
#' This is only reported if each list in \code{...} was generated with \code{pval.type="any"} in \code{\link{combineMarkers}}.
#' \item \code{p.value}, the largest p-value of each gene across all DataFrames for that group.
#' This is replaced by \code{log.p.value} if p-values in \code{...} are log-transformed.
#' \item \code{FDR}, the BH-adjusted value of \code{p.value}.
#' This is replaced by \code{log.FDR} if p-values in \code{...} are log-transformed.
#' }
#' 
#' @details
#' The combined statistics are designed to favor a gene that is highly ranked in each of the individual test results.
#' This is highly conservative and aims to identify robust DE that is significant under all testing schemes.
#'
#' A combined \code{Top} value of T indicates that the gene is among the top T genes of one or more pairwise comparisons
#' in each of the testing schemes.
#' (We can be even more aggressive if the individual results were generated with a larger \code{min.prop} value.)
#' In effect, a gene can only achieve a low \code{Top} value if it is consistently highly ranked in each test.
#' If \code{sorted=TRUE}, this is used to order the genes in the output DataFrame.
#'
#' The combined \code{p.value} is effectively the result of applying an intersection-union test to the per-test results.
#' This will only be low if the gene has a low p-value in each of the test results.
#' If \code{sorted=TRUE} and \code{Top} is not present, this will be used to order the genes in the output DataFrame.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{findMarkers}} and \code{\link{combineMarkers}}, to generate elements in \code{...}.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay, only using k-means for convenience.
#' kout <- kmeans(t(logcounts(sce)), centers=4) 
#'
#' tout <- findMarkers(sce, groups=kout$cluster, direction="up")
#' wout <- findMarkers(sce, groups=kout$cluster, direction="up", test="wilcox")
#'
#' combined <- multiMarkerStats(t=tout, wilcox=wout)
#' colnames(combined[[1]])
#' 
#' @export
#' @importFrom S4Vectors DataFrame SimpleList I
multiMarkerStats <- function(..., repeated=NULL, sorted=TRUE) {
    all.methods <- list(...)
    nmethods <- length(all.methods)
    if (is.null(names(all.methods))) {
        stop("elements in '...' must be named")
    }

    # Checking groups within each element.
    all.methods <- .check_element_names(all.methods, 
        FUN=names, SUB=function(x, i) x[i], 
        msg="names inside each list")
    group.names <- names(all.methods[[1]])

    # Iterating over each group and combining the statistics.
    output <- vector("list", length(group.names))
    names(output) <- group.names

    for (i in seq_along(group.names)) {
        collected <- vector("list", nmethods)
        for (j in seq_len(nmethods)) {
            collected[[j]] <- all.methods[[j]][[i]]
        }

        collected <- .check_element_names(collected, 
            FUN=rownames, SUB=function(x, i) x[i,,drop=FALSE],
            msg=sprintf("row names for '%s'", group.names[i]))
        gene.names <- rownames(collected[[1]])

        # Pulling out the redundant columns before counting the number of columns.
        redundancies <- list()
        for (r in repeated) {
            curval <- NULL
            for (j in seq_len(nmethods)) {
                if (!is.null(collected[[j]][[r]])) {
                    curval <- collected[[j]][[r]]
                    collected[[j]][[r]] <- NULL
                }
            }
            redundancies[[r]] <- curval
        }

        ncols <- vapply(collected, ncol, 0L)
        if (length(unique(ncols))!=1L) {
            stop(sprintf("different numbers of columns for '%s'", group.names[i]))
        }
        ncols <- ncols[1]

        collated <- DataFrame(row.names=gene.names)
        if (length(redundancies)) {
            redundancies <- do.call(DataFrame, redundancies)
            collated <- cbind(collated, redundancies)
        }

        # Interleave relevant statistics. We assume they're in the same order across DFs.
        interleaved <- vector("list", nmethods)
        for (j in seq_len(nmethods)) {
            current <- lapply(collected[[j]], I)
            names(current) <- paste0(names(all.methods)[j], ".", names(current))
            interleaved[[j]] <- current
        }
    
        flip <- as.integer(matrix(seq_len(nmethods*ncols), nrow=nmethods, byrow=TRUE))
        interleaved <- unlist(interleaved, recursive=FALSE)[flip]
        interleaved <- do.call(DataFrame, interleaved)

        # Computing combined statistics.
        all.top <- .extract_columns(collected, "Top")
        all.p <- .extract_columns(collected, "p.value")
        all.lp <- .extract_columns(collected, "log.p.value")

        if (!is.null(all.top)) {
            collated$Top <- do.call(pmax, all.top)
        }
        if (!is.null(all.p)) {
            stats <-collated$p.value <- do.call(pmax, all.p)
            collated$FDR <- p.adjust(collated$p.value, method="BH")
        } else if (!is.null(all.lp)) {
            stats <- collated$log.p.value <- do.call(pmax, all.lp)
            collated$log.FDR <- .logBH(collated$log.p.value)
        } else {
            stop("p-value field should be 'p.value' or 'log.p.value'")
        }

        # Assembling the final result.
        collated <- cbind(collated, interleaved)
        if (sorted) {
            o <- order(if (!is.null(collated$Top)) collated$Top else stats)
            collated <- collated[o,,drop=FALSE]
        }
        output[[i]] <- collated
    }

    SimpleList(output)
}

.check_element_names <- function(things, FUN, SUB, msg) {
    .names <- FUN(things[[1]])
    s.names <- sort(.names)
    for (i in seq_along(things)) {
        cur.names <- FUN(things[[i]])
        if (is.null(cur.names)) {
            stop(paste(msg, "should be non-NULL"))
        }
        if (!identical(s.names, sort(FUN(things[[i]])))) {
            stop(paste(msg, "should be the same"))
        }
        things[[i]] <- SUB(things[[i]], .names)
    }

    things
}

.extract_columns <- function(things, col) {
    for (j in seq_along(things)) {
        current <- things[[j]]
        things[j] <- list(if (col %in% colnames(current)) current[,col] else NULL)
    }

    decision <- unique(vapply(things, is.null, TRUE))
    if (length(decision)!=1L) {
        stop(sprintf("either all or no elements should have '%s'", col))
    }

    if (decision) NULL else things
}
