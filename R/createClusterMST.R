#' Minimum spanning trees on cluster centroids
#'
#' Perform basic trajectory analyses with minimum spanning trees (MST) computed on cluster centroids,
#' based on the methodology in the \pkg{TSCAN} package.
#'
#' @param centers A numeric matrix of cluster centroids where each \emph{row} represents a cluster 
#' and each column represents a dimension (usually a PC or another low-dimensional embedding).
#' Each row should be named with the cluster name.
#' @param mst A \link{graph} object containing a MST, typically the output of \code{createClusterMST(centers)}.
#' For \code{connectClusterMSTNodes}, the MST may be computed from a different \code{centers}.
#' @param combined Logical scalar indicating whether a single data.frame of edge coordinates should be returned.
#' @param x A numeric matrix of per-cell coordinates where each \emph{row} represents a cell
#' and each column represents a dimension (again, usually a low-dimensional embedding).
#' @param ids A character vector of length equal to the number of cells,
#' specifying the cluster to which each cell is assigned.
#' @param start A string specifying the starting node from which to compute pseudotimes.
#' Defaults to an arbitrarily chosen node of degree 1.
#'
#' @details
#' These functions represent some basic utilities for a simple trajectory analysis 
#' based on the algorithm in the \pkg{TSCAN} package.
#'
#' \code{createClusterMST} builds a MST where each node is a cluster centroid and 
#' each edge is weighted by the Euclidean distance between centroids.
#' This represents the most parsimonious explanation for a particular trajectory
#' and has the advantage of being directly intepretable with respect to any pre-existing clusters.
#' 
#' \code{connectClusterMST} provides the coordinates of the start and end of every edge.
#' This is mostly useful for plotting purposes in \code{\link{segments}} or the equivalent \pkg{ggplot2} functionality.
#' We suggest using \code{\link{aggregateAcrossCells}} to obtain \code{centers} for multiple low-dimensional results at once.
#'
#' \code{orderClusterMST} will map each cell to the closest edge involving the cluster to which it is assigned.
#' (Here, edges are segments terminated by their nodes, so some cells may simply be mapped to the edge terminus.)
#' It will then calculate the distance of that cell along the MST from the starting node specified by \code{start}.
#' This distance represents the pseudotime for that cell and can be used in further quantitative analyses.
#'
#' @section Breaking down the pseudotime matrix:
#' The pseudotimes are returned as a matrix where each row corresponds to cell in \code{x} 
#' and each column corresponds to a path through the MST from \code{start} to all nodes of degree 1.
#' (If \code{start} is itself a node of degree 1, then paths are only considered to all other such nodes.)
#' This format is inspired by that from the \pkg{slingshot} package and provides a compact representation of branching events.
#'
#' Each branching event in the MST results in a new path and thus a new column in the pseudotime matrix.
#' For any given row in this matrix, entries are either \code{NA} or they are identical.
#' This reflects the fact that multiple paths will share a section of the MST for which the pseudotimes are the same.
#' 
#' The starting node in \code{start} is \emph{completely arbitrarily chosen} by \code{orderClusterMST},
#' as directionality is impossible to infer from the expression matrix alone.
#' However, it is often possible to use prior biological knowledge to pick an appropriate cluster as the starting node.
#'
#' @return 
#' \code{createClusterMST} returns a \link{graph} object containing an MST computed on \code{centers}.
#'
#' \code{connectClusterMST} returns, by default, a data.frame containing the start and end coordinates of segments representing all the edges in \code{mst}.
#' If \code{combined=FALSE}, a list of two data.frames is returned where corresponding rows represent the start and end coordinates of the same edge.
#'
#' \code{orderClusterMST} returns a numeric matrix containing the pseudotimes of all cells (rows) across all paths (columns) through \code{mst}.
#'
#' @author Aaron Lun
#' @references
#' Ji Z and Ji H (2016).
#' TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
#' \emph{Nucleic Acids Res.} 44, e117
#'
#' @examples
#' # Mocking up a Y-shaped trajectory.
#' centers <- rbind(c(0,0), c(0, -1), c(1, 1), c(-1, 1))
#' rownames(centers) <- seq_len(nrow(centers))
#' clusters <- sample(nrow(centers), 1000, replace=TRUE)
#' cells <- centers[clusters,]
#' cells <- cells + rnorm(length(cells), sd=0.5)
#'
#' # Creating the MST first:
#' mst <- createClusterMST(centers)
#' plot(mst)
#'
#' # Also plotting the MST on top of existing visualizations:
#' edges <- connectClusterMST(centers, mst, combined=FALSE)
#' plot(cells[,1], cells[,2], col=clusters)
#' segments(edges$start$dim1, edges$start$dim2, edges$end$dim1, 
#'      edges$end$dim2, lwd=5)
#' 
#' # Finally, obtaining pseudo-time orderings.
#' ordering <- orderClusterMST(cells, clusters, centers, mst)
#' unified <- rowMeans(ordering, na.rm=TRUE)
#' plot(cells[,1], cells[,2], col=topo.colors(21)[cut(unified, 21)], pch=16)
#' 
#' @name createClusterMST
NULL

#' @export
#' @rdname createClusterMST
#' @importFrom igraph graph.adjacency minimum.spanning.tree
createClusterMST <- function(centers) {
    dmat <- dist(centers)
    dmat <- as.matrix(dmat)
    g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
    minimum.spanning.tree(g)
}

#' @export
#' @rdname createClusterMST
#' @importFrom Matrix which
#' @importFrom igraph V
connectClusterMST <- function(centers, mst, combined=TRUE) {
    pairs <- which(mst[] > 0, arr.ind=TRUE)
    pairs <- pairs[pairs[,1] > pairs[,2],,drop=FALSE]

    vnames <- names(V(mst))
    group <- paste0(vnames[pairs[,1]], "--", vnames[pairs[,2]])

    if (is.null(colnames(centers))) {
        colnames(centers) <- sprintf("dim%i", seq_len(ncol(centers)))
    }

    L <- vnames[pairs[,1]]
    R <- vnames[pairs[,2]]
    L <- data.frame(edge=group, centers[L,,drop=FALSE])
    R <- data.frame(edge=group, centers[R,,drop=FALSE])
    rownames(L) <- rownames(R) <- NULL

    if (combined) {
        rbind(L, R)
    } else {
        list(start=L, end=R)
    }
}

#' @export
#' @rdname createClusterMST
#' @importFrom igraph V degree adjacent_vertices
orderClusterMST <- function(x, ids, centers, mst, start=NULL) {
    if (is.null(start)) {
        start <- names(V(mst)[degree(mst)==1])[1]
    }
    if (!all(ids %in% rownames(centers))) {
        stop("all 'ids' must be in 'rownames(centers)'")
    }

    collated <- list()
    latest <- as.character(start)
    parents <- NA_character_ 
    progress <- list(rep(NA_real_, length(ids)))
    cumulative <- 0

    while (length(latest)) {
        new.latest <- new.parents <- character(0)
        new.progress <- list()
        new.cumulative <- numeric(0)

        for (i in seq_along(latest)) {
            curnode <- latest[i]
            all.neighbors <- names(adjacent_vertices(mst, curnode, mode="all")[[1]])
            cur.ids <- ids==curnode 

            mapped <- .map2edges(x[cur.ids,,drop=FALSE], center=centers[curnode,], 
                edge.ends=centers[all.neighbors,,drop=FALSE], previous=parents[i])
            edge.len <- mapped$dist
            pseudo <- mapped$pseudo

            cum.dist <- cumulative[i] + edge.len

            collected.progress <- list()
            for (j in seq_along(pseudo)) {
                sofar <- progress[[i]] # yes, the 'i' here is deliberate.
                sofar[cur.ids] <- pseudo[[j]] + cum.dist
                collected.progress[[j]] <- sofar
            }

            all.children <- setdiff(all.neighbors, parents[i])
            if (length(all.children)==0) {
                collated[[curnode]] <- collected.progress[[1]]
            } else {
                new.latest <- c(new.latest, all.children)
                new.parents <- c(new.parents, rep(curnode, length(all.children)))
                new.progress <- c(new.progress, collected.progress)
                new.cumulative <- c(new.cumulative, rep(cum.dist, length(all.children)))
            }
        }

        latest <- new.latest
        parents <- new.parents
        progress <- new.progress
        cumulative <- new.cumulative
    }
    
    do.call(cbind, collated)
}

.map2edges <- function(points, center, edge.ends, previous) {
    all.distances <- list()
    all.pseudo <- list()
    edge.len <- list()

    centered <- t(t(points) - center)

    # Computing distance of each point from each edge.
    # Edges defined from 'center' to 'edge.ends'.
    for (i in rownames(edge.ends)) {
        edge.end <- edge.ends[i,]
        delta <- edge.end - center
        max.d <- sqrt(sum(delta^2))
        delta <- delta/max.d

        proj <- as.numeric(centered %*% delta)
        proj <- pmax(0, pmin(proj, max.d))
        mapped <- outer(proj, delta)

        dist <- sqrt(rowSums((centered - mapped)^2))
        all.distances[[i]] <- dist
        all.pseudo[[i]] <- proj
        edge.len[[i]] <- max.d
    }

    all.distances <- do.call(cbind, all.distances)
    all.pseudo <- do.call(cbind, all.pseudo)
    chosen <- colnames(all.distances)[max.col(-all.distances, ties.method="first")]

    # Flipping the distance of points to the previous node,
    # in order to enforce a directional pseudotime.
    dist.previous <- 0
    if (!is.na(previous)) {
        on.previous <- chosen==previous
        dist.previous <- edge.len[[previous]]
        previous.proj <- -all.pseudo[on.previous,previous,drop=FALSE]

        if (all(on.previous)) {
            return(list(dist=dist.previous, pseudo=list(previous.proj)))
        }
    }

    # Filling out the branches, where points are NA for a branch's
    # pseudotime if they were assigned to another branch.
    output <- list()
    for (leftover in setdiff(rownames(edge.ends), previous)) {
        empty <- rep(NA_real_, nrow(points))
        if (!is.na(previous)) {
            empty[on.previous] <- previous.proj
        }
        current <- chosen==leftover
        empty[current] <- all.pseudo[current,leftover]
        output[[leftover]] <- empty
    }

    list(dist=dist.previous, pseudo=output)
}
