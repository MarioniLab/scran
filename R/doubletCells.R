#' @importFrom scater librarySizeFactors normalizeCounts
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom Matrix rowMeans
#' @importFrom stats median
#' @importFrom BiocNeighbors findKNN findNeighbors queryNeighbors queryKNN
.doublet_cells <- function(x, size.factors.norm=NULL, size.factors.content=NULL,
    k=50, subset.row=NULL, niters=max(10000, ncol(x)), block=10000, 
    d=50, approximate=FALSE, irlba.args=list(), 
    force.match=FALSE, force.k=20, force.ndist=3,
    BNPARAM=NULL, BPPARAM=SerialParam())
# Simulates doublets and uses a mutual nearest-neighbour approach to match them to real cells.
#
# written by Aaron Lun
# created 18 July 2018
{
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (is.null(size.factors.norm)) {
        size.factors.norm <- librarySizeFactors(x)
    }
    if (!is.null(size.factors.content)) {
        x <- normalizeCounts(x, size.factors.content, return_log=FALSE, centre_size_factors=FALSE)
        size.factors.norm <- size.factors.norm/size.factors.content
    }

    y <- normalizeCounts(x, size.factors.norm, centre_size_factors=FALSE)

    # Running the SVD.
    svd.out <- .centered_SVD(t(y), max.rank=d, approximate=approximate, extra.args=irlba.args, 
        keep.left=TRUE, keep.right=TRUE)
    pcs <- .svd_to_pca(svd.out, ncomp=d, named=FALSE)
    sim.pcs <- .spawn_doublet_pcs(x, size.factors.norm, V=svd.out$v, centers=rowMeans(y), niters=niters, block=block)

    # Force doublets to nearest neighbours in the original data set.
    pre.pcs <- precluster(pcs)
    if (force.match) {
        closest <- queryKNN(query=sim.pcs, k=force.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, precomputed=pre.pcs)
        sim.pcs <- .compute_tricube_average(pcs, closest$index, closest$distance, ndist=force.ndist)
    }

    # Computing densities, using a distance computed from the kth nearest neighbor.
    self.dist <- findKNN(precomputed=pre.pcs, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.index=FALSE)$distance
    dist2nth <- pmax(1e-8, median(self.dist[,ncol(self.dist)]))

    self.dist <- findNeighbors(precomputed=pre.pcs, threshold=dist2nth, BPPARAM=BPPARAM, get.index=FALSE)$distance
    sim.dist <- queryNeighbors(sim.pcs, query=pcs, threshold=dist2nth, BPPARAM=BPPARAM, get.index=FALSE)$distance

    rel.dens <- bpmapply(FUN=function(self, sim, limit) {
        sum((1 - (sim/limit)^3)^3)/sum((1 - (self/limit)^3)^3)^2
    }, self=self.dist, sim=sim.dist, limit=dist2nth, BPPARAM=BPPARAM)

    rel.dens/(niters/ncol(x))
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
        sim.sf <- size.factors[left] + size.factors[right]
        sim.y <- normalizeCounts(sim.x, sim.sf, centre_size_factors=FALSE)

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
setGeneric("doubletCells", function(x, ...) standardGeneric("doubletCells"))

#' @export
setMethod("doubletCells", "ANY", .doublet_cells)

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
setMethod("doubletCells", "SingleCellExperiment", function(x, size.factors.norm=NA, ..., subset.row=NULL, assay.type="counts", get.spikes=FALSE) {
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    if (length(size.factors.norm)==1L && is.na(size.factors.norm)) {
        size.factors.norm <- sizeFactors(x)
    }
    .doublet_cells(assay(x, i=assay.type), size.factors.norm=size.factors.norm, ..., subset.row=subset.row)
})
