#' @importFrom scater librarySizeFactors normalize
#' @importFrom BiocGenerics "sizeFactors<-" sizeFactors
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom Matrix crossprod rowMeans
#' @importFrom stats median
#' @importFrom kmknn findKNN findNeighbors queryNeighbors
.doublet_cells <- function(x, d=50, k=20, approximate=FALSE, subset.row=NULL, multiple=1, block=10000, BPPARAM=SerialParam())
# Simulates doublets and uses a mutual nearest-neighbour approach to match them to real cells.
# 
# written by Aaron Lun
# created 18 July 2018 
{
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    sce <- SingleCellExperiment(list(counts=x))
    sizeFactors(sce) <- librarySizeFactors(x)
    sce <- normalize(sce, return_log=TRUE)
    y <- logcounts(sce, withDimnames=FALSE)

    # Running the SVD.
    args <- list(y=t(y), max.rank=d, value="lowrank")
    if (approximate) {
        svdfun <- .irlba_svd
        args <- c(args, irlba.args)
    } else {
        svdfun <- .full_svd
    }
    svd.out <- do.call(svdfun, args)
    pcs <- sweep(svd.out$u, 2L, svd.out$d, FUN="*", check.margin = FALSE)

    # Simulating doublets.
    collected <- list()
    counter <- 1L
    current <- 0L
    niters <- round(multiple * ncol(x))

    while (current < niters) {
        to.make <- min(block, niters - current)
        sim <- x[,sample(ncol(x), to.make, replace=TRUE),drop=FALSE] + x[,sample(ncol(x), to.make, replace=TRUE),drop=FALSE]
            
		# Normalizing by their library sizes.
        sim.sce <- SingleCellExperiment(list(counts=sim))
        sizeFactors(sim.sce) <- librarySizeFactors(sim.sce)
        sim.sce <- normalize(sim.sce, return_log=TRUE)
        sim.y <- logcounts(sim.sce, withDimnames=FALSE)
       
        # Projecting onto the PC space of the original data.
        sim.pcs <- crossprod(sim.y, svd.out$v)
        collected[[counter]] <- sim.pcs
        counter <- counter + 1L
        current <- current + block
    }
    
    sim.pcs <- do.call(rbind, collected)
    all.means <- rowMeans(y)
    correction <- colSums(all.means * svd.out$v)
    sim.pcs <- sweep(sim.pcs, 2L, correction, FUN="-", check.margin=FALSE)
    
    # Computing densities, using a distance computed from the kth nearest neighbor.
    pre.pcs <- precluster(pcs)
    self.dist <- findKNN(precomputed=pre.pcs, k=k, BPPARAM=BPPARAM, get.index=FALSE)$distance
    dist2nth <- median(self.dist[,ncol(self.dist)])

    self.dist <- findNeighbors(precomputed=pre.pcs, threshold=dist2nth, BPPARAM=BPPARAM, get.index=FALSE)$distance
    sim.dist <- queryNeighbors(sim.pcs, query=pcs, threshold=dist2nth, BPPARAM=BPPARAM, get.index=FALSE)$distance
    rel.dens <- bpmapply(FUN=function(self, sim, limit) {
        sum((1 - (sim/limit)^3)^3)/sum((1 - (self/limit)^3)^3)
    }, self=self.dist, sim=sim.dist, limit=dist2nth, BPPARAM=BPPARAM)

    rel.dens/multiple
}
