#' Detect doublet cells
#' 
#' Identify potential doublet cells based on simulations of putative doublet expression profiles.
#' 
#' @param x A numeric matrix-like object of count values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#' 
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param size.factors.norm A numeric vector of size factors for normalization of \code{x} prior to PCA and distance calculations.
#' If \code{NULL}, defaults to size factors derived from the library sizes of \code{x}.
#' 
#' For the SingleCellExperiment method, the default values are taken from \code{\link{sizeFactors}(x)}, if they are available.
#' @param size.factors.content A numeric vector of size factors for RNA content normalization of \code{x} prior to simulating doublets.
#' This is orthogonal to the values in \code{size.factors.norm}, see Details.
#' @param k An integer scalar specifying the number of nearest neighbours to use to determine the bandwidth for density calculations.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param niters An integer scalar specifying how many simulated doublets should be generated.
#' @param block An integer scalar controlling the rate of doublet generation, to keep memory usage low.
#' @param d An integer scalar specifying the number of components to retain after the PCA.
#' @param force.match A logical scalar indicating whether remapping of simulated doublets to original cells should be performed.
#' @param force.k An integer scalar specifying the number of neighbours to use for remapping if \code{force.match=TRUE}.
#' @param force.ndist A numeric scalar specifying the bandwidth for remapping if \code{force.match=TRUE}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' This should be an algorithm supported by \code{\link{findNeighbors}}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, if \code{d} is not \code{NA}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the neighbour searches should be parallelized.
#' @param ... For the generic, additional arguments to pass to specific methods.
#' 
#' For the SingleCellExperiment method, additional arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values contain the count matrix. 
#' @param get.spikes See \code{?"\link{scran-gene-selection}"}.
#' 
#' @return 
#' A numeric vector of doublet scores for each cell in \code{x}.
#' 
#' @details
#' This function simulates doublets by adding the count vectors for two randomly chosen cells in \code{x}.
#' For each original cell, we compute the density of neighboring simulated doublets and compare it to the density of neighboring original cells.
#' Genuine doublets should have a high density of simulated doublets relative to the density of its neighbourhood.
#' Thus, the doublet score for each cell is defined as the ratio of densities of simulated doublets to the (squared) density of the original cells.
#' 
#' Densities are calculated in low-dimensional space after a PCA on the log-normalized expression matrix of \code{x}.
#' Simulated doublets are projected into the low-dimensional space using the rotation vectors computed from the original cells.
#' A tricube kernel is used to compute the density around each cell.
#' The bandwidth of the kernel is set to the median distance to the \code{k} nearest neighbour across all cells.
#' 
#' The two size factor arguments have different roles:
#' \itemize{
#' \item \code{size.factors.norm} contains the size factors to be used for normalization prior to PCA and distance calculations.
#' This defaults to the values returned by \code{\link{librarySizeFactors}} but can be explicitly set to ensure that the low-dimensional space is consistent with that in the rest of the analysis.
#' \item \code{size.factors.content} is much more important, and represents the size factors that preserve RNA content differences.
#' This is usually computed from spike-in RNA and ensures that the simulated doublets have the correct ratio of contributions from the original cells.
#' }
#' It is possible to set both of these arguments as they are orthogonal to each other.
#' Setting \code{size.factors.content} will not affect the calculation of log-normalized expression values from \code{x}.
#' Conversely, setting \code{size.factors.norm} will not affect the ratio in which cells are added together when simulating doublets.
#' 
#' If \code{force.match=TRUE}, simulated doublets will be remapped to the nearest neighbours in the original data.
#' This is done by taking the (tricube-weighted) average of the PC scores for the \code{force.k} nearest neighbors.
#' The tricube bandwidth for remapping is chosen by taking the median distance and multiplying it by \code{force.ndist}, to protect against later neighbours that might be outliers.
#' The aim is to adjust for unknown differences in RNA content that would cause the simulated doublets to be systematically displaced from their true locations.
#' However, it may also result in spuriously high scores for single cells that happen to be close to a cluster of simulated doublets.
#' @author
#' Aaron Lun
#' 
#' @examples
#' # Mocking up an example.
#' ngenes <- 1000
#' mu1 <- 2^rnorm(ngenes)
#' mu2 <- 2^rnorm(ngenes)
#' mu3 <- 2^rnorm(ngenes)
#' mu4 <- 2^rnorm(ngenes)
#' 
#' counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes) # Pure type 1
#' counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes) # Pure type 2
#' counts.3 <- matrix(rpois(ngenes*100, mu3), nrow=ngenes) # Pure type 3
#' counts.4 <- matrix(rpois(ngenes*100, mu4), nrow=ngenes) # Pure type 4
#' counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes) # Doublets (1 & 2)
#' 
#' counts <- cbind(counts.1, counts.2, counts.3, counts.4, counts.m)
#' clusters <- rep(1:5, c(rep(100, 4), ncol(counts.m)))
#' 
#' # Find potential doublets.
#' scores <- doubletCells(counts)
#' boxplot(split(log10(scores), clusters))
#' 
#' @references
#' Lun ATL (2018).
#' Detecting doublet cells with \emph{scran}.
#' \url{https://ltla.github.io/SingleCellThoughts/software/doublet_detection/bycell.html}
#'
#' @seealso
#' \code{\link{doubletCluster}}, to detect doublet clusters.
#' @name doubletCells
NULL

#' @importFrom scater librarySizeFactors normalizeCounts
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom BiocParallel SerialParam bpmapply bpstart bpisup bpstop
#' @importFrom Matrix rowMeans
#' @importFrom stats median
#' @importFrom BiocNeighbors findKNN findNeighbors queryNeighbors queryKNN buildIndex
#' @importFrom methods is
.doublet_cells <- function(x, size.factors.norm=NULL, size.factors.content=NULL,
    k=50, subset.row=NULL, niters=max(10000, ncol(x)), block=10000, 
    d=50, force.match=FALSE, force.k=20, force.ndist=3,
    BNPARAM=KmknnParam(), BSPARAM=ExactParam(), BPPARAM=SerialParam())
{
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (is.null(size.factors.norm)) {
        size.factors.norm <- librarySizeFactors(x)
    }

    # Manually controlling the size factor centering here to ensure the final counts are on the same scale.
    size.factors.norm <- size.factors.norm/mean(size.factors.norm)
    if (!is.null(size.factors.content)) {
        x <- normalizeCounts(x, size.factors.content, log=FALSE, center_size_factors=FALSE)
        size.factors.norm <- size.factors.norm/size.factors.content
    }
    y <- normalizeCounts(x, size.factors.norm, center_size_factors=FALSE)

    if (!bpisup(BPPARAM)){ 
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Running the SVD.
    svd.out <- .centered_SVD(t(y), max.rank=d, keep.left=TRUE, keep.right=TRUE, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    pcs <- .svd_to_pca(svd.out, ncomp=d, named=FALSE)
    sim.pcs <- .spawn_doublet_pcs(x, size.factors.norm, V=svd.out$v, centers=rowMeans(y), niters=niters, block=block)

    # Force doublets to nearest neighbours in the original data set.
    pre.pcs <- buildIndex(pcs, BNPARAM=BNPARAM)
    if (force.match) {
        closest <- queryKNN(query=sim.pcs, k=force.k, BNINDEX=pre.pcs, BPPARAM=BPPARAM)
        sim.pcs <- .compute_tricube_average(pcs, closest$index, closest$distance, ndist=force.ndist)
    }

    # Computing densities, using a distance computed from the kth nearest neighbor.
    self.dist <- findKNN(BNINDEX=pre.pcs, k=k, BPPARAM=BPPARAM, get.distance=FALSE, get.index=FALSE)
    dist2nth <- pmax(1e-8, median(self.dist))

    self.dist <- findNeighbors(threshold=dist2nth, BNINDEX=pre.pcs, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.index=FALSE)$distance
    self.prop <- lengths(self.dist)/ncol(x)
    sim.dist <- queryNeighbors(sim.pcs, query=pcs, threshold=dist2nth, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.index=FALSE)$distance
    sim.prop <- lengths(sim.dist)/niters

    sim.prop / self.prop^2
}

#' @importFrom Matrix crossprod
#' @importFrom scater normalizeCounts
.spawn_doublet_pcs <- function(x, size.factors, V, centers, niters=10000L, block=10000L) {
    collected <- list()
    counter <- 1L
    current <- 0L
    mean.correction <- colSums(centers * V)

    while (current < niters) {
        to.make <- min(block, niters - current)
        left <- sample(ncol(x), to.make, replace=TRUE)
        right <- sample(ncol(x), to.make, replace=TRUE)
        sim.x <- x[,left,drop=FALSE] + x[,right,drop=FALSE]

        # Do not center, otherwise the simulated doublets will always have higher normalized counts
        # than actual doublets (as the latter will have been normalized to the level of singlets).
        sim.sf <- size.factors[left] + size.factors[right]
        sim.y <- normalizeCounts(sim.x, sim.sf, center_size_factors=FALSE)

        # Projecting onto the PC space of the original data.
        sim.pcs <- crossprod(sim.y, V)
        sim.pcs <- sweep(sim.pcs, 2L, mean.correction, FUN="-", check.margin=FALSE)
        collected[[counter]] <- sim.pcs
        counter <- counter + 1L
        current <- current + block
    }

    do.call(rbind, collected)
}

##############################
# S4 method definitions here #
##############################

#' @export
#' @rdname doubletCells
setGeneric("doubletCells", function(x, ...) standardGeneric("doubletCells"))

#' @export
#' @rdname doubletCells
setMethod("doubletCells", "ANY", .doublet_cells)

#' @export
#' @rdname doubletCells
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
setMethod("doubletCells", "SingleCellExperiment", function(x, size.factors.norm=sizeFactors(x), ..., 
    subset.row=NULL, assay.type="counts", get.spikes=FALSE) 
{
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .doublet_cells(assay(x, i=assay.type), size.factors.norm=size.factors.norm, ..., subset.row=subset.row)
})
