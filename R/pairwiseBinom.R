#' Perform pairwise binomial tests
#'
#' Perform pairwise binomial tests between groups of cells, 
#' possibly after blocking on uninteresting factors of variation.
#' 
#' @param x A numeric matrix-like object of counts,
#' where each column corresponds to a cell and each row corresponds to a gene.
#' @param clusters A vector of cluster identities for all cells.
#' @param block A factor specifying the blocking level for each cell.
#' @param direction A string specifying the direction of effects to be considered for each cluster.
#' @param log.p A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param gene.names A character vector of gene names with one value for each row of \code{x}.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param threshold Numeric scalar specifying the value below which a gene is presumed to be not expressed.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how parallelization should be performed across genes.
#' 
#' @details
#' This function performs exact binomial tests to identify marker genes between pairs of clusters.
#' Here, the null hypothesis is that the proportion of cells expressing a gene is the same between clusters.
#' A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of clusters.
#' This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.
#' 
#' Effect sizes for each comparison are reported as log2-fold changes 
#' in the proportion of expressing cells in one cluster over the proportion in another cluster.
#' Large log-FCs correspond to large relative differences in these proportions,
#' where the sign indicates the direction of the change in proportions.
#' We add a pseudo-count that squeezes the log-FCs towards zero, to avoid undefined values when one proportion is zero.
#' 
#' \code{x} can be a count matrix or any transformed counterpart where zeroes remain zero and non-zeroes remain non-zero.
#' This is true of any scaling normalization and monotonic transformation like the log-transform.
#' If the transformation breaks this rule, some adjustment of \code{threshold} is necessary.
#' 
#' A consequence of the transformation-agnostic behaviour of this function is that it will not respond to normalization.
#' Differences in library size will not be considered by this function.
#' However, this is not necessarily problematic for marker gene detection -
#' users can treat this as \emph{retaining} information about the total RNA content, analogous to spike-in normalization.
#' 
#' @section Direction and magnitude of the effect:
#' If \code{direction="any"}, two-sided binomial tests will be performed for each pairwise comparisons between clusters.
#' For other \code{direction}, one-sided tests in the specified direction will be used to compute p-values for each gene.
#' This can be used to focus on genes that are upregulated in each cluster of interest, which is often easier to interpret.
#' 
#' In practice, the two-sided test is approximated by combining two one-sided tests using a Bonferroni correction.
#' This is done for various logistical purposes;
#' it is also the only way to combine p-values across blocks in a direction-aware manner.
#' As a result, the two-sided p-value reported here will not be the same as that from \code{\link{binom.test}}.
#' In practice, they are usually similar enough that this is not a major concern.
#' 
#' To interpret the setting of \code{direction}, consider the DataFrame for cluster X, in which we are comparing to another cluster Y.
#' If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in cluster X compared to Y.
#' If \code{direction="down"}, genes will only be significant if they are downregulated in cluster X compared to Y.
#' See \code{?\link{binom.test}} for more details on the interpretation of one-sided Wilcoxon rank sum tests.
#'
#' The magnitude of the log-fold change in the proportion of expressing cells can also be tested by setting \code{lfc}.
#' By default, \code{lfc=0} meaning that we will reject the null upon detecting any difference in proportions.
#' If this is set to some other positive value, the null hypothesis will change depending on \code{direction}:
#' \itemize{
#' \item If \code{direction="any"}, the null hypothesis is that the true log-fold change in proportions is either \code{-lfc} or \code{lfc} with equal probability.
#' A two-sided p-value is computed against this composite null.
#' \item If \code{direction="up"}, the null hypothesis is that the true log-fold change is \code{lfc}, and a one-sided p-value is computed.
#' \item If \code{direction="down"}, the null hypothesis is that the true log-fold change is \code{-lfc}, and a one-sided p-value is computed.
#' }
#' 
#' @section Blocking on uninteresting factors:
#' If \code{block} is specified, binomial tests are performed between clusters within each level of \code{block}.
#' For each pair of clusters, the p-values for each gene across 
#' all levels of \code{block} are combined using Stouffer's weighted Z-score method.
#' 
#' The weight for the p-value in a particular level of \code{block} is defined as \eqn{N_x + N_y},
#' where \eqn{N_x} and \eqn{N_y} are the number of cells in clusters X and Y, respectively, for that level. 
#' This means that p-values from blocks with more cells will have a greater contribution to the combined p-value for each gene.
#' 
#' When combining across batches, one-sided p-values in the same direction are combined first.
#' Then, if \code{direction="any"}, the two combined p-values from both directions are combined.
#' This ensures that a gene only receives a low overall p-value if it changes in the same direction across batches.
#'
#' When comparing two clusters, blocking levels are ignored if no p-value was reported, e.g., if there were insufficient cells for a cluster in a particular level. 
#' If all levels are ignored in this manner, the entire comparison will only contain \code{NA} p-values and a warning will be emitted.
#' 
#' @return
#' A list is returned containing \code{statistics} and \code{pairs}.
#' 
#' The \code{statistics} element is itself a list of \linkS4class{DataFrame}s.
#' Each DataFrame contains the statistics for a comparison between a pair of clusters,
#' including the overlap proportions, p-values and false discovery rates.
#' 
#' The \code{pairs} element is a DataFrame with one row corresponding to each entry of \code{statistics}.
#' This contains the fields \code{first} and \code{second}, 
#' specifying the two clusters under comparison in the corresponding DataFrame in \code{statistics}.
#' 
#' In each DataFrame in \code{statistics}, the log-fold change represents the log-ratio of the proportion of expressing cells in the \code{first} cluster compared to the expressing proportion in the \code{second} cluster.
#' 
#' @author
#' Aaron Lun
#' 
#' @references
#' Whitlock MC (2005). 
#' Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
#' \emph{J. Evol. Biol.} 18, 5:1368-73.
#' 
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(example.sce)), centers=2) 
#' 
#' # Vanilla application:
#' out <- pairwiseBinom(logcounts(example.sce), clusters=kout$cluster)
#' out
#' 
#' # Directional and with a minimum log-fold change:
#' out <- pairwiseBinom(logcounts(example.sce), clusters=kout$cluster, 
#'     direction="up", lfc=1)
#' out
#'
#' @seealso
#' \code{\link{binom.test}}, on which this function is based.
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
pairwiseBinom <- function(x, clusters, block=NULL, direction=c("any", "up", "down"),
    log.p=FALSE, gene.names=rownames(x), subset.row=NULL, threshold=1e-8, lfc=0, 
    BPPARAM=SerialParam())
{
    ncells <- ncol(x)
    clusters <- as.factor(clusters)
    if (length(clusters)!=ncells) {
        stop("length of 'clusters' does not equal 'ncol(x)'")
    }
    if (nlevels(clusters) < 2L) {
        stop("need at least two unique levels in 'clusters'")
    }

    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    if (is.null(gene.names)) {
        gene.names <- subset.row
    } else if (length(gene.names)!=nrow(x)) {
        stop("length of 'gene.names' is not equal to the number of rows")
    } else {
        gene.names <- gene.names[subset.row]
    }

    direction <- match.arg(direction)
    results <- .blocked_binom(x, subset.row, clusters, block=block, direction=direction, 
        gene.names=gene.names, log.p=log.p, threshold=threshold, lfc=lfc, BPPARAM=BPPARAM)

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names), use.names=FALSE)
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply SerialParam bpisup bpstart bpstop
#' @importFrom stats pbinom
.blocked_binom <- function(x, subset.row, clusters, block=NULL, direction="any", gene.names=NULL, log.p=TRUE, 
	threshold=1e-8, lfc=0, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# binomial tests between pairs of clusters within each blocking level.
{
    if (is.null(block)) {
        block <- list(`1`=seq_len(ncol(x)))
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- split(seq_along(block), block)
    }

    # Choosing the parallelization strategy.
    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row, wout)

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Computing across blocks.
    clust.vals <- levels(clusters)
    nblocks <- length(block)
    all.nzero <- all.n <- vector("list", nblocks)

    for (b in seq_along(block)) {
        chosen <- block[[b]]
        cur.clusters <- clusters[chosen]
        all.n[[b]] <- as.vector(table(cur.clusters))
        names(all.n[[b]]) <- clust.vals
        
        cur.groups <- split(chosen, cur.clusters)
        raw.nzero <- bplapply(by.core, FUN=.compute_nzero_stat, x=x, by.group=cur.groups, 
            threshold=threshold, BPPARAM=BPPARAM)
        cons.nzero <- do.call(rbind, raw.nzero)
        colnames(cons.nzero) <- clust.vals
        all.nzero[[b]] <- cons.nzero
    }

    if (lfc==0) {
        STATFUN <- .generate_nolfc_binom(all.n, all.nzero)
    } else {
        STATFUN <- .generate_lfc_binom(all.n, all.nzero, direction, lfc)
    }

    .pairwise_blocked_template(x, clust.vals, nblocks=length(block), direction=direction, 
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, effect.name="logFC")
}

##########################
### Internal functions ###
##########################

#' @importFrom stats pbinom
.generate_nolfc_binom <- function(all.n, all.nzero) {
    force(all.n)
    force(all.nzero)

    function(b, host, target) {    
        host.nzero <- all.nzero[[b]][,host]
        target.nzero <- all.nzero[[b]][,target]
        host.n <- all.n[[b]][[host]]
        target.n <- all.n[[b]][[target]]

        p <- host.n/(host.n + target.n)
        size <- host.nzero + target.nzero
        effect <- .compute_binom_effect(host.nzero, host.n, target.nzero, target.n)

        # Not exactly equal to binom.test(); the two-sided p-value from binom.test()
        # cannot be performed by any combination of the one-sided p-values. This 
        # makes it impossible to behave with directional Stouffer's method in 
        # .pairwise_blocking_template(), and just generally gums up the works.
        list(
            forward=effect, 
            reverse=-effect,
            weight=as.double(host.n) + as.double(target.n),
            valid=host.n > 0L && target.n > 0L,
            left=pbinom(host.nzero, size, p, log.p=TRUE),
            right=pbinom(host.nzero - 1, size, p, lower.tail=FALSE, log.p=TRUE)
        )
    }
}

# Log-fold change in proportions, mimic edgeR::cpm().
.compute_binom_effect <- function(host.nzero, host.n, target.nzero, target.n) {
    mean.lib <- mean(c(host.n, target.n))
    pseudo.host <- 1 * host.n/mean.lib
    pseudo.target <- 1 * target.n/mean.lib
    unname(log2((host.nzero + pseudo.host)/(host.n + 2 * pseudo.host))
        - log2((target.nzero + pseudo.target)/(target.n + 2 * pseudo.target)))
}

#' @importFrom stats pbinom
.generate_lfc_binom <- function(all.n, all.nzero, direction, lfc) {
    force(all.n)
    force(all.nzero)
    force(direction)
    fold <- 2^lfc

    function(b, host, target) {    
        host.nzero <- all.nzero[[b]][,host]
        target.nzero <- all.nzero[[b]][,target]
        host.n <- all.n[[b]][[host]]
        target.n <- all.n[[b]][[target]]

        # Log-fold change in the ratios.
        p.left <- host.n/fold / (target.n + host.n/fold)
        p.right <- host.n*fold / (target.n + host.n*fold)
        size <- host.nzero + target.nzero
        effect <- .compute_binom_effect(host.nzero, host.n, target.nzero, target.n)

        output <- list(
            forward=effect, 
            reverse=-effect,
            weight=as.double(host.n) + as.double(target.n),
            valid=host.n > 0L && target.n > 0L
        )

        left.lower <- pbinom(host.nzero, size, p.left, log.p=TRUE)
        right.upper <- pbinom(host.nzero - 1, size, p.right, lower.tail=FALSE, log.p=TRUE)

        if (direction=="any") {
            left.upper <- pbinom(host.nzero, size, p.right, log.p=TRUE)
            right.lower <- pbinom(host.nzero - 1, size, p.left, lower.tail=FALSE, log.p=TRUE)

            # Here, the null hypothesis is that the shift is evenly distributed at 50%
            # probability for -lfc and lfc, hence we take the average of the two p-values.
            output$left <- .add_log_values(left.lower, left.upper) - log(2)
            output$right <- .add_log_values(right.lower, right.upper) - log(2)
        } else {
            output$left <- left.lower
            output$right <- right.upper
        }

        output
    }
}

#' @importFrom scater nexprs
.compute_nzero_stat <- function(x, by.group, rows, threshold) {
    collected <- lapply(by.group, FUN=function(s) {
        nexprs(x, subset_col=s, subset_row=rows, byrow=TRUE, detection_limit=threshold)
    })
    do.call(cbind, collected)
}
