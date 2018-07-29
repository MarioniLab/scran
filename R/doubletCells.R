#' @importFrom scater librarySizeFactors normalizeMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom Matrix rowMeans
#' @importFrom stats median
#' @importFrom kmknn findKNN findNeighbors queryNeighbors
.doublet_cells <- function(x, size.factors.norm=NULL, size.factors.content=NULL,
    k=50, subset.row=NULL, niters=max(10000, ncol(x)), block=10000, force.match=FALSE, force.k=20,
    d=50, approximate=FALSE, irlba.args=list(), BPPARAM=SerialParam())
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
        x <- normalizeMatrix(x, size.factors.content, return_log=FALSE, centre_size_factors=FALSE)
        size.factors.norm <- size.factors.norm/size.factors.content
    }

    y <- normalizeMatrix(x, size.factors.norm, centre_size_factors=FALSE)

    # Running the SVD.
    svd.out <- .centered_SVD(t(y), max.rank=d, approximate=approximate, extra.args=irlba.args, 
        keep.left=TRUE, keep.right=TRUE)
    pcs <- .svd_to_pca(svd.out, ncomp=d, named=FALSE)
    sim.pcs <- .spawn_doublet_pcs(x, size.factors.norm, V=svd.out$v, centers=rowMeans(y), niters=niters, block=block)

    # Computing densities, using a distance computed from the kth nearest neighbor.
    pre.pcs <- precluster(pcs)
    self.dist <- findKNN(precomputed=pre.pcs, k=k, BPPARAM=BPPARAM, get.index=FALSE)$distance
    dist2nth <- median(self.dist[,ncol(self.dist)])

    if (force.match) {
        sim.pcs <- .tricube_weighted_remapping(sim.pcs, pcs, precomputed=pre.pcs, k=force.k, BPPARAM=BPPARAM)        
    }

    self.dist <- findNeighbors(precomputed=pre.pcs, threshold=dist2nth, BPPARAM=BPPARAM, get.index=FALSE)$distance
    sim.dist <- queryNeighbors(sim.pcs, query=pcs, threshold=dist2nth, BPPARAM=BPPARAM, get.index=FALSE)$distance
    rel.dens <- bpmapply(FUN=function(self, sim, limit) {
        sum((1 - (sim/limit)^3)^3)/sum((1 - (self/limit)^3)^3)
    }, self=self.dist, sim=sim.dist, limit=dist2nth, BPPARAM=BPPARAM)

    rel.dens/(niters/ncol(x))
}

#' @importFrom Matrix crossprod
#' @importFrom scater normalizeMatrix
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
        sim.y <- normalizeMatrix(sim.x, sim.sf, centre_size_factors=FALSE)

        # Projecting onto the PC space of the original data.
        sim.pcs <- crossprod(sim.y, V)
        sim.pcs <- sweep(sim.pcs, 2L, mean.correction, FUN="-", check.margin=FALSE)
        collected[[counter]] <- sim.pcs
        counter <- counter + 1L
        current <- current + block
    }

    do.call(rbind, collected)
}

#' @importFrom kmknn queryKNN
.tricube_weighted_remapping <- function(simulated, original, ...)
# Computing tricube-weighted coordinates for remapping simulated cells to their closest originals.
{
    closest <- queryKNN(query=simulated, X=original, ...)
    max.dist <- closest$distance[,ncol(closest$distance)]
    rel.dist <- closest$distance / max.dist

    tricube <- (1 - rel.dist^3)^3
    weight <- tricube/rowSums(tricube)
    output <- matrix(0, nrow(simulated), ncol(simulated))
    for (kdx in seq_len(ncol(closest$index))) {
        output <- output + original[closest$index[,kdx],,drop=FALSE] * weight[,kdx]
    }
    
    return(output)
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
