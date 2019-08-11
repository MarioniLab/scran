#' Perform pairwise t-tests
#' 
#' Perform pairwise Welch t-tests between groups of cells, possibly after blocking on uninteresting factors of variation.
#' 
#' @param x A numeric matrix-like object of normalized log-expression values, 
#' where each column corresponds to a cell and each row corresponds to an endogenous gene.
#' @param clusters A vector of cluster identities for all cells.
#' @param block A factor specifying the blocking level for each cell.
#' @param design A numeric matrix containing blocking terms for uninteresting factors.
#' Note that these should not be confounded with \code{clusters} or contain an intercept, see Details.
#' @param direction A string specifying the direction of log-fold changes to be considered for each cluster.
#' @param lfc A positive numeric scalar specifying the log-fold change threshold to be tested against.
#' @param std.lfc A logical scalar indicating whether log-fold changes should be standardized.
#' @param log.p A logical scalar indicating if log-transformed p-values/FDRs should be returned.
#' @param gene.names A character vector of gene names with one value for each row of \code{x}.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes.
#' 
#' @details
#' This function performs Welch t-tests to identify differentially expressed genes (DEGs) between pairs of clusters.
#' The aim is to use the DEGs to determine cluster identity based on expression of marker genes with known biological activity.
#' A list of tables is returned where each table contains the statistics for all genes for a comparison between each pair of clusters.
#' This can be examined directly or used as input to \code{\link{combineMarkers}} for marker gene detection.
#' 
#' The Welch t-test is simple, fast and performs reasonably well for single-cell count data (Soneson and Robinson, 2018).
#' However, if one of the clusters contains fewer than two cells, no p-value will be reported for comparisons involving that cluster.
#' A warning will also be raised about insufficient degrees of freedom (d.f.) in such cases.
#' 
#' When \code{log.p=TRUE}, the log-transformed p-values and FDRs are reported using the natural base.
#' This is useful in cases with many cells such that reporting the p-values directly would lead to double-precision underflow.
#' 
#' @section Direction and magnitude of the log-fold change:
#' Log-fold changes are reported as differences in the values of \code{x}.
#' Thus, all log-fold changes have the same base as whatever was used to perform the log-transformation in \code{x}.
#' If \code{\link[scater]{normalize}} was used, this would be base 2.
#' 
#' If \code{direction="any"}, two-sided tests will be performed for each pairwise comparisons between clusters.
#' Otherwise, one-sided tests in the specified direction will be used to compute p-values for each gene.
#' This can be used to focus on genes that are upregulated in each cluster of interest, which is often easier to interpret.
#' 
#' To interpret the setting of \code{direction}, consider the DataFrame for cluster X, in which we are comparing to another cluster Y.
#' If \code{direction="up"}, genes will only be significant in this DataFrame if they are upregulated in cluster X compared to Y.
#' If \code{direction="down"}, genes will only be significant if they are downregulated in cluster X compared to Y.
#' 
#' The magnitude of the log-fold changes can also be tested by setting \code{lfc}.
#' By default, \code{lfc=0} meaning that we will reject the null upon detecting any differential expression.
#' If this is set to some other positive value, the null hypothesis will change depending on \code{direction}:
#' \itemize{
#' \item If \code{direction="any"}, the null hypothesis is that the true log-fold change is either \code{-lfc} or \code{lfc} with equal probability.
#' A two-sided p-value is computed against this composite null.
#' \item If \code{direction="up"}, the null hypothesis is that the true log-fold change is \code{lfc}, and a one-sided p-value is computed.
#' \item If \code{direction="down"}, the null hypothesis is that the true log-fold change is \code{-lfc}, and a one-sided p-value is computed.
#' }
#' This is similar to the approach used in \code{\link[limma:eBayes]{treat}} and allows users to focus on genes with strong log-fold changes.
#' 
#' If \code{std.lfc=TRUE}, the log-fold change for each gene is standardized by the variance.
#' When the Welch t-test is being used, this is equivalent to Cohen's d.
#' Standardized log-fold changes may be more appealing for visualization as it avoids large fold changes due to large variance.
#' The choice of \code{std.lfc} does not affect the calculation of the p-values.
#' 
#' @section Blocking on uninteresting factors:
#' If \code{block} is specified, t-tests are performed between clusters within each level of \code{block}.
#' For each pair of clusters, the p-values for each gene across all levels of \code{block} are combined using Stouffer's weighted Z-score method.
#' The reported log-fold change for each gene is also a weighted average of log-fold changes across levels.
#' 
#' The weight for a particular level is defined as \eqn{(1/N_x + 1/N_y)^{-1}}, 
#' where \eqn{Nx} and \eqn{Ny} are the number of cells in clusters X and Y, respectively, for that level. 
#' This is inversely proportional to the expected variance of the log-fold change, provided that all clusters and blocking levels have the same variance.
#' 
#' % In theory, a better weighting scheme would be to use the estimated standard error of the log-fold change to compute the weight.
#' % This would be more responsive to differences in variance between blocking levels, focusing on levels with low variance and high power.
#' % However, this is not safe in practice as genes with many zeroes can have very low standard errors, dominating the results inappropriately.
#' 
#' When comparing two clusters, blocking levels are ignored if no p-value was reported, e.g., if there were insufficient cells for a cluster in a particular level. 
#' This includes levels that contain fewer than two cells for either cluster, as this cannot yield a p-value from the Welch t-test.
#' If all levels are ignored in this manner, the entire comparison will only contain \code{NA} p-values and a warning will be emitted.
#' 
#' @section Regressing out unwanted factors:
#' If \code{design} is specified, a linear model is instead fitted to the expression profile for each gene.
#' This linear model will include the \code{clusters} as well as any blocking factors in \code{design}.
#' A t-test is then performed to identify DEGs between pairs of clusters, using the values of the relevant coefficients and the gene-wise residual variance.
#' Note that \code{design} must be full rank when combined with the \code{clusters} terms, i.e., there should not be any confounding variables.
#' Similarly, any intercept column should be removed beforehand.
#' 
#' We recommend using \code{block} instead of \code{design} for uninteresting categorical factors of variation.
#' The former accommodates differences in the variance of expression in each cluster via Welch's t-test.
#' As a result, it is more robust to misspecification of the clusters, as misspecified clusters (and inflated variances) do not affect the inferences for other clusters.
#' Use of \code{block} also avoids assuming additivity of effects between the blocking factors and the cluster identities.
#' 
#' Nonetheless, use of \code{design} is unavoidable when blocking on real-valued covariates.
#' It is also useful for ensuring that log-fold changes/p-values are computed for comparisons between all pairs of clusters
#' (assuming that \code{design} is not confounded with the cluster identities).
#' This may not be the case with \code{block} if a pair of clusters never co-occur in a single blocking level. 
#' 
#' @return
#' A list is returned containing \code{statistics} and \code{pairs}.
#' 
#' The \code{statistics} element is itself a list of \linkS4class{DataFrame}s.
#' Each DataFrame contains the statistics for a comparison between a pair of clusters,
#' including the log-fold changes, p-values and false discovery rates.
#' 
#' The \code{pairs} element is a DataFrame where each row corresponds to an entry of \code{statistics}.
#' This contains the fields \code{first} and \code{second}, 
#' specifying the two clusters under comparison in the corresponding DataFrame in \code{statistics}.
#' 
#' In each DataFrame in \code{statistics}, the log-fold change represents the change in the \code{first} cluster compared to the \code{second} cluster.
#' 
#' @author
#' Aaron Lun
#' 
#' @references
#' Whitlock MC (2005). 
#' Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. 
#' \emph{J. Evol. Biol.} 18, 5:1368-73.
#' 
#' Soneson C and Robinson MD (2018). 
#' Bias, robustness and scalability in single-cell differential expression analysis. 
#' \emph{Nat. Methods}
#' 
#' Lun ATL (2018).
#' Comments on marker detection in \emph{scran}.
#' \url{https://ltla.github.io/SingleCellThoughts/software/marker_detection/comments.html}
#'
#' @seealso
#' \code{\link{t.test}}, on which this function is based.
#' 
#' @examples
#' data(example.sce)
#'
#' # Any clustering method is okay.
#' kout <- kmeans(t(logcounts(example.sce)), centers=3) 
#' 
#' # Vanilla application:
#' out <- pairwiseTTests(logcounts(example.sce), clusters=kout$cluster)
#' out
#' 
#' # Directional with log-fold change threshold:
#' out <- pairwiseTTests(logcounts(example.sce), clusters=kout$cluster, 
#'     direction="up", lfc=0.2)
#' out
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
pairwiseTTests <- function(x, clusters, block=NULL, design=NULL, direction=c("any", "up", "down"),
    lfc=0, std.lfc=FALSE, log.p=FALSE, gene.names=rownames(x), subset.row=NULL, BPPARAM=SerialParam())
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

    if (!is.null(block) && !is.null(design)) {
        stop("cannot specify both 'block' and 'design'")
    } else if (!is.null(design)) {
        results <- .fit_lm_internal(x, subset.row, clusters, design=design, direction=direction, lfc=lfc, 
            std.lfc=std.lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
    } else {
        results <- .test_block_internal(x, subset.row, clusters, block=block, direction=direction, lfc=lfc, 
            std.lfc=std.lfc, gene.names=gene.names, log.p=log.p, BPPARAM=BPPARAM)
    }

    first <- rep(names(results), lengths(results))
    second <- unlist(lapply(results, names), use.names=FALSE)
    results <- unlist(results, recursive=FALSE, use.names=FALSE)
    names(results) <- NULL
    list(statistics=results, pairs=DataFrame(first=first, second=second))
}

###########################################################
# Internal functions (blocking)
###########################################################

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply SerialParam
.test_block_internal <- function(x, subset.row, clusters, block=NULL, direction="any", lfc=0, std.lfc=FALSE,
    gene.names=NULL, log.p=TRUE, BPPARAM=SerialParam())
# This looks at every level of the blocking factor and performs
# t-tests between pairs of clusters within each blocking level.
{
    nclusters <- nlevels(clusters)
    if (is.null(block)) {
        all.clusters <- factor(as.integer(clusters), seq_len(nclusters))
        nblocks <- 1L
    } else {
        if (length(block)!=ncol(x)) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- factor(block)
        nblocks <- nlevels(block)
        all.clusters <- factor(as.integer(clusters) + (as.integer(block) - 1L) * nclusters, seq_len(nblocks*nclusters))
    }

    # Calculating the statistics for each block.
    all.blocks <- split(seq_along(all.clusters) - 1L, all.clusters)
    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row, wout)
    by.core <- .split_matrix_by_workers(x, by.core)

    raw.stats <- bplapply(by.core, FUN=compute_blocked_stats_none, bygroup=all.blocks, BPPARAM=BPPARAM)
    all.means <- do.call(rbind, lapply(raw.stats, FUN=function(x) t(x[[1]])))
    all.vars <- do.call(rbind, lapply(raw.stats, FUN=function(x) t(x[[2]])))
    all.n <- table(all.clusters)

    clust.vals <- levels(clusters)
    out.s2 <- out.means <- out.n <- vector("list", nblocks)
    for (b in seq_len(nblocks)) {
        chosen <- (b-1L) * nclusters + seq_len(nclusters)
        means <- all.means[,chosen,drop=FALSE]
        sigma2 <- all.vars[,chosen,drop=FALSE]
        cur.n <- as.vector(all.n[chosen])

        colnames(means) <- colnames(sigma2) <- names(cur.n) <- clust.vals
        out.means[[b]] <- means
        out.s2[[b]] <- sigma2
        out.n[[b]] <- cur.n 
    }

    # Running through all pairs of comparisons.
    STATFUN <- function(b, host, target) {
        host.n <- out.n[[b]][host]
        target.n <- out.n[[b]][target]
        host.s2 <- out.s2[[b]][,host]
        target.s2 <- out.s2[[b]][,target]
        t.out <- .get_t_test_stats(host.s2=host.s2, target.s2=target.s2, host.n=host.n, target.n=target.n)

        cur.err <- t.out$err
        cur.df <- t.out$test.df
        cur.lfc <- out.means[[b]][,host] - out.means[[b]][,target]
        p.out <- .run_t_test(cur.lfc, cur.err, cur.df, thresh.lfc=lfc, direction=direction)

        effect.size <- cur.lfc
        if (std.lfc) {
            # Computing Cohen's D.
            pooled.s2 <- ((host.n - 1) * host.s2 + (target.n - 1) * target.s2)/(target.n + host.n - 2)
            effect.size <- effect.size / sqrt(pooled.s2)
        }
        
        list(forward=effect.size, reverse=-effect.size, left=p.out$left, right=p.out$right, 

            # Weights are inversely proportional to the squared error of the log-fold change,
            # _assuming equal variance_ across blocks and groups for simplicity.
            # (Using the actual variance is dangerous as some blocks have zero variance
            # with a defined log-fold change, if there are many zeroes.)
            weight=1/(1/host.n + 1/target.n),

            # Indicating that there is no valid test statistic if the df is zero.
            valid=all(!is.na(cur.df))
        )
    }

    .pairwise_blocked_template(x, clust.vals, nblocks, direction=direction,
        gene.names=gene.names, log.p=log.p, STATFUN=STATFUN, effect.name="logFC")
}

.get_t_test_stats <- function(host.s2, target.s2, host.n, target.n)
# Computes the error variance of the log-fold change and the
# degrees of freedom for the t-distribution.
{
    host.df <- max(0L, host.n - 1L)
    target.df <- max(0L, target.n - 1L)

    # Avoid unlikely but potential problems with discreteness.
    host.s2 <- pmax(host.s2, 1e-8)
    target.s2 <- pmax(target.s2, 1e-8)

    if (host.df > 0L && target.df > 0L) {
        # Perform Welch's t-test here.
        host.err <- host.s2/host.n
        target.err <- target.s2/target.n
        cur.err <- host.err + target.err
        cur.df <- cur.err^2 / (host.err^2/host.df + target.err^2/target.df)
    } else {
        cur.err <- cur.df <- NA_real_
    }
    return(list(err=cur.err, test.df=cur.df))
}

###########################################################
# Internal functions (linear modelling)
###########################################################

#' @importFrom stats model.matrix
#' @importFrom limma lmFit contrasts.fit
#' @importFrom BiocParallel bplapply SerialParam
.fit_lm_internal <- function(x, subset.row, clusters, design, direction="any", lfc=0, std.lfc=FALSE,
    gene.names=NULL, log.p=TRUE, BPPARAM=SerialParam())
# This fits a linear model to each gene and performs a t-test for
# differential expression between clusters.
{
    full.design <- model.matrix(~0 + clusters)
    clust.vals <- levels(clusters)
    colnames(full.design) <- clust.vals

    if (nrow(design)!=ncol(x)) {
        stop("'nrow(design)' is not equal to 'ncol(x)'")
    }
    full.design <- cbind(full.design, design) # Other linear dependencies will trigger errors in .ranksafe_QR.

    # Getting coefficient estimates (need to undo column pivoting to get the actual estimates).
    QR <- .ranksafe_qr(full.design)
    resid.df <- nrow(full.design) - ncol(full.design)
    if (resid.df <= 0L) {
        stop("no residual d.f. in design matrix for variance estimation")
    }

    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row-1L, wout)
    raw.stats <- bplapply(by.core, FUN=fit_linear_model, qr=QR$qr, qraux=QR$qraux, exprs=x, get_coefs=TRUE, BPPARAM=BPPARAM)
    coefficients <- do.call(cbind, lapply(raw.stats, "[[", i=1))
    coefficients[QR$pivot,] <- coefficients
    sigma2 <- unlist(lapply(raw.stats, "[[", i=3))
    sigma2 <- pmax(sigma2, 1e-8) # avoid unlikely but possible problems with discreteness.

    # Running through every pair of clusters.
    out.stats <- .create_output_container(clust.vals)
    ngenes <- length(subset.row)
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)

    for (h in seq_along(clust.vals)) {
        host <- clust.vals[h]

        # Computing standard errors (via limma, to avoid having to manually calculate standard errors).
        con <- matrix(0, ncol(full.design), h-1L)
        diag(con) <- -1
        con[h,] <- 1
        lfit2 <- contrasts.fit(lfit, con)

        # Computing log-fold changes, t-statistics and p-values on a per-contrast basis.
        # This _could_ be vectorised, but it's less confusing to do it like this,
        # and there's not much speed gain to be had from vectorizing over contrasts.
        ref.coef <- coefficients[h,]
        for (tdex in seq_len(h-1L)) {
            target <- clust.vals[tdex]
            cur.lfc <- ref.coef - coefficients[tdex,]

            test.out <- .run_t_test(cur.lfc, lfit2$stdev.unscaled[tdex]^2*sigma2, resid.df, thresh.lfc=lfc, direction=direction)
            hvt.p <- .choose_leftright_pvalues(test.out$left, test.out$right, direction=direction)
            tvh.p <- .choose_leftright_pvalues(test.out$right, test.out$left, direction=direction)

            if (std.lfc) {
                # Computing Cohen's D.
                cur.lfc <- cur.lfc / sqrt(sigma2)
            }
            out.stats[[host]][[target]] <- .create_full_stats(cur.lfc, p=hvt.p, gene.names=gene.names, log.p=log.p)
            out.stats[[target]][[host]] <- .create_full_stats(-cur.lfc, p=tvh.p, gene.names=gene.names, log.p=log.p)
        }
    }

    out.stats
}

###########################################################
# Internal functions (t-test calculations)
###########################################################

#' @importFrom stats pt
.run_t_test <- function(cur.lfc, cur.err, cur.df, thresh.lfc=0, direction="any")
# This runs the t-test given the relevant statistics, regardless of how
# they were computed (i.e., within blocks, or in a linear model).
{
    thresh.lfc <- abs(thresh.lfc)
    if (thresh.lfc==0) {
        cur.t <- cur.lfc/sqrt(cur.err)
        left <- pt(cur.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
        right <- pt(cur.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
    } else {
        upper.t <- (cur.lfc - thresh.lfc)/sqrt(cur.err)
        lower.t <- (cur.lfc + thresh.lfc)/sqrt(cur.err)

        left.lower <- pt(lower.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
        right.upper <- pt(upper.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)

        if (direction=="any") {
            # Using the TREAT method, which tests against a null where the log-fold change is
            # takes values at the extremes of [-thresh, thresh]. The null probability is 50%
            # distributed across both extremes, hence the -log(2) at the end.
            left.upper <- pt(upper.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
            right.lower <- pt(lower.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
            left <- .add_log_values(left.upper, left.lower) - log(2)
            right <- .add_log_values(right.upper, right.lower) - log(2)
        } else {
            # For one-sided tests, testing against the lfc threshold in the specified direction.
            # Note that if direction='up', only 'right' is used; nonetheless, 'left' is still calculated
            # using 'lower.t' so as to allow quick calculation of the p-value for the reversed contrast,
            # by simply swapping 'left' and 'right'. The same applies when direction='down'.
            left <- left.lower
            right <- right.upper
        }
    }
    return(list(left=left, right=right))
}

.add_log_values <- function(x, y)
# This performs log(exp(x) + exp(y)) in a numerically
# stable fashion, i.e., without actually running the 'exp's.
{
    pmax(x, y) + log1p(exp(-abs(x-y)))
}
