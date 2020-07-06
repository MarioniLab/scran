#' Minimum spanning trees on cluster centroids
#'
#' Perform basic trajectory analyses with minimum spanning trees (MST) computed on cluster centroids,
#' based on the methodology in the \pkg{TSCAN} package.
#'
#' @param centers A numeric matrix of cluster centroids where each \emph{row} represents a cluster 
#' and each column represents a dimension (usually a PC or another low-dimensional embedding).
#' Each row should be named with the cluster name.
#' @param outgroup A logical scalar indicating whether an outgroup should be inserted to split unrelated trajectories.
#' Alternatively, a numeric scalar specifying the distance threshold to use for this splitting.
#' @param outscale A numeric scalar specifying the scaling to apply to the median distance between centroids
#' to define the threshold for outgroup splitting.
#' Only used if \code{outgroup=TRUE}.
#' @param mst A \link{graph} object containing a MST, typically the output of \code{createClusterMST(centers)}.
#' For \code{connectClusterMSTNodes}, the MST may be computed from a different \code{centers}.
#' @param combined Logical scalar indicating whether a single data.frame of edge coordinates should be returned.
#' @param x A numeric matrix of per-cell coordinates where each \emph{row} represents a cell
#' and each column represents a dimension (again, usually a low-dimensional embedding).
#' @param ids A character vector of length equal to the number of cells,
#' specifying the cluster to which each cell is assigned.
#' @param start A character vector specifying the starting node from which to compute pseudotimes in each component of \code{mst}.
#' Defaults to an arbitrarily chosen node of degree 1 or lower in each component.
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
#' @section Introducing an outgroup:
#' If \code{outgroup=TRUE}, we add an outgroup to avoid constructing a trajectory between \dQuote{unrelated} clusters.
#' This is done by adding an extra row/column to the distance matrix corresponding to an artificial outgroup cluster,
#' where the distance to all of the other real clusters is set to \eqn{\omega/2}.
#' Large jumps in the MST between real clusters that are more distant than \eqn{\omega} will then be rerouted through the outgroup,
#' allowing us to break up the MST into multiple subcomponents by removing the outgroup.
#'
#' The default \eqn{\omega} value is computed by constructing the MST from the original distance matrix,
#' computing the median edge length in that MST, and then scaling it by \code{outscale}.
#' This adapts to the magnitude of the distances and the internal structure of the dataset
#' while also providing some margin for variation across cluster pairs.
#' Alternatively, \code{outgroup} can be set to a numeric scalar in which case it is used directly as \eqn{\omega}.
#'
#' @section Confidence on the edges:
#' For the MST, we obtain a measure of the confidence in each edge by computing the distance gained if that edge were not present.
#' Ambiguous parts of the tree will be less penalized from deletion of an edge, manifesting as a small distance gain.
#' In contrast, parts of the tree with clear structure will receive a large distance gain upon deletion of an obvious edge.
#'
#' For each edge, we divide the distance gain by the length of the edge to normalize for cluster resolution.
#' This avoids overly penalizing edges in parts of the tree involving broad clusters
#' while still retaining sensitivity to detect distance gain in overclustered regions.
#' As an example, a normalized gain of unity for a particular edge means that its removal
#' requires an alternative path that increases the distance travelled by that edge's length.
#'
#' The normalized gain is reported as the \code{"gain"} attribute in the edges of the MST from \code{\link{createClusterMST}}.
#' Note that the \code{"weight"} attribute represents the edge length.
#' 
#' @section Interpreting the pseudotime matrix:
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
#' @seealso
#' \code{\link{quickPseudotime}}, a wrapper to quickly perform these calculations.
#' @name createClusterMST
NULL

#' @export
#' @rdname createClusterMST
#' @importFrom igraph graph.adjacency minimum.spanning.tree delete_vertices E
#' @importFrom stats median dist
createClusterMST <- function(centers, outgroup=FALSE, outscale=3) {
    dmat <- dist(centers)
    dmat <- as.matrix(dmat)

    if (!isFALSE(outgroup)) {
        if (!is.numeric(outgroup)) {
            g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
            mst <- minimum.spanning.tree(g)
            med <- median(E(mst)$weight)
            outgroup <- med * outscale
        }

        old.d <- rownames(dmat)

        # Divide by 2 so rerouted distance between cluster pairs is 'outgroup'.
        dmat <- rbind(cbind(dmat, outgroup/2), outgroup/2) 
        diag(dmat) <- 0

        special.name <- strrep("x", max(nchar(old.d))+1L)
        rownames(dmat) <- colnames(dmat) <- c(old.d, special.name)
    }

    g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
    mst <- minimum.spanning.tree(g)
    mst <- .estimate_edge_confidence(mst, g)

    if (!isFALSE(outgroup)) {
        mst <- delete_vertices(mst, special.name)
    }

    mst 
}

#' @importFrom igraph minimum.spanning.tree E E<- ends get.edge.ids delete.edges
.estimate_edge_confidence <- function(mst, g) {
    edges <- E(mst)
    ends <- ends(mst, edges)
    reweight <- numeric(length(edges))

    for (i in seq_along(edges)) {
        id <- get.edge.ids(g, ends[i,])        
        g.copy <- delete.edges(g, id)
        mst.copy <- minimum.spanning.tree(g.copy)
        reweight[i] <- sum(E(mst.copy)$weight)
    }

    W <- edges$weight
    total <- sum(W)
    offset <- min(W)
    E(mst)$gain <- (reweight - total)/(W + offset/1e8)
    mst
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
#' @importFrom igraph V degree adjacent_vertices components
orderClusterMST <- function(x, ids, centers, mst, start=NULL) {
    comp <- components(mst)$membership
    by.comp <- split(names(comp), comp)
    if (is.null(start)) {
        candidates <- names(V(mst)[degree(mst) <= 1])
        start <- vapply(by.comp, function(b) intersect(b, candidates)[1], "")
    } else {
        start <- as.character(start)
        for (b in by.comp) {
            if (length(intersect(b, start))!=1) {
                stop("'start' must have one cluster in each component of 'mst'")
            }
        }
    }

    if (!all(ids %in% rownames(centers))) {
        stop("all 'ids' must be in 'rownames(centers)'")
    }

    collated <- list()
    latest <- start
    nstarts <- length(latest)
    parents <- rep(NA_character_, nstarts)
    progress <- rep(list(rep(NA_real_, length(ids))), nstarts)
    cumulative <- numeric(nstarts)

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
    if (nrow(edge.ends)==0L) {
        # Return _something_ with a 1-cluster input.
        return(list(dist=0, pseudo=numeric(nrow(points))))
    }

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
    best <- max.col(-all.distances, ties.method="first")
    chosen <- colnames(all.distances)[best]

    # Flipping the distance of points to the previous node,
    # in order to enforce a directional pseudotime.
    dist.previous <- 0
    if (!is.na(previous)) {
        on.previous <- chosen==previous 
        dist.previous <- edge.len[[previous]]
        previous.proj <- -all.pseudo[on.previous,previous,drop=FALSE]

        if (ncol(all.distances)==1) {
            return(list(dist=dist.previous, pseudo=list(previous.proj)))
        }
    }

    # If any distances are zero, the corresponding cells are considered to be
    # shared with all paths, as they are assigned right at the branch point.
    dist <- all.pseudo[cbind(seq_along(best), best)]
    in.everyone <- dist==0

    # Filling out the branches, where points are NA for a branch's
    # pseudotime if they were assigned to another branch.
    output <- list()
    for (leftover in setdiff(rownames(edge.ends), previous)) {
        empty <- rep(NA_real_, nrow(points))
        empty[in.everyone] <- 0
        if (!is.na(previous)) {
            empty[on.previous] <- previous.proj
        }
        current <- chosen==leftover 
        empty[current] <- all.pseudo[current,leftover]
        output[[leftover]] <- empty
    }

    list(dist=dist.previous, pseudo=output)
}
