# This tests the findMarkers function.
# require(scran); require(testthat); source("test-markers.R")

REFFUN <- function(y, grouping, pval.type="any", direction="any", full.stats=FALSE) 
# A reference function using the t.test function.
{ 
    output <- findMarkers(y, grouping, pval.type=pval.type, direction=direction, full.stats=full.stats)
    grouping <- factor(grouping)
    clust.vals <- levels(grouping)
    alt.hyp <- switch(direction, any="two.sided", up="greater", down="less")

    for (host in clust.vals) {
        collected.lfc <- collected.p <- list()
        host.y <- y[,grouping==host,drop=FALSE]

        for (target in setdiff(clust.vals, host)) {
            target.y <- y[,grouping==target,drop=FALSE]

            if (ncol(host.y)>1L && ncol(target.y)>1L) {
                pval <- numeric(nrow(y))
                for (i in seq_along(pval)) {
                    pval[i] <- t.test(host.y[i,], target.y[i,], alternative=alt.hyp, var.equal=FALSE)$p.value
                }
            } else {
                pval <- rep(NA_real_, nrow(y))
            }
            
            collected.lfc[[paste0("logFC.", target)]] <- rowMeans(host.y) - rowMeans(target.y)
            collected.p[[target]] <- pval
        }
        
        # Compiling the requested ordering. 
        combined.p <- do.call(cbind, collected.p)
        if (pval.type=="any") { 
            all.ranks <- lapply(collected.p, rank, ties.method="first", na.last="keep")
            reordered <- order(do.call(pmin, c(all.ranks, na.rm=TRUE)), 
                               do.call(pmin, c(collected.p, na.rm=TRUE)))
            pval <- apply(combined.p, 1, FUN=function(x) { min(p.adjust(x, method="BH"), na.rm=TRUE) })

            # Checking that the ordering is consistent with the ranking for each p.value set.
            standard <- output[[host]]
            orderings <- lapply(collected.p, FUN=function(p) { order(p, na.last=NA) })
            for (best in seq(1, 101, by=5)) { 
                best.in.each <- lapply(orderings, "[", i=seq_len(best))
                expect_identical(sort(rownames(standard)[standard$Top <= best]),
                                 sort(rownames(y)[unique(unlist(best.in.each))]))
            }
        } else {
            pval <- apply(combined.p, 1, FUN=max, na.rm=TRUE)
            reordered <- order(pval)
        }
       
        # Running the requested checks.
        cur.out <- output[[host]]
        expect_identical(rownames(cur.out), rownames(y)[reordered])
        adj.p <- p.adjust(pval, method="BH")
        expect_equal(adj.p[reordered], cur.out$FDR)

        if (!full.stats) {
            expect_equal(as.matrix(cur.out[,-(1:2)]), do.call(cbind, collected.lfc)[reordered,,drop=FALSE])
        } else {
            for (target in setdiff(clust.vals, host)) {
                cur.stats <- cur.out[,paste0("stats.", target)]
                expect_equal(cur.stats$logFC, unname(collected.lfc[[paste0("logFC.", target)]][reordered]))
                expect_equal(cur.stats$p.value, collected.p[[target]][reordered])
                expect_equal(cur.stats$FDR, p.adjust(collected.p[[target]][reordered], method="BH"))
            }
        }
    }  
    return(TRUE)
}

set.seed(70000)
ncells <- 200
ngenes <- 250
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

set.seed(7000001)
test_that("findmarkers works as expected without blocking or design matrices", {
    clust <- kmeans(t(exprs(X)), centers=3)
    clusters <- as.factor(clust$cluster)

    # Vanilla runs.
    REFFUN(exprs(X), clusters)
    REFFUN(exprs(X), clusters, direction="up")
    REFFUN(exprs(X), clusters, direction="down")
    REFFUN(exprs(X), clusters, full.stats=TRUE)

    # Checking that the IUT calculation is correct.
    REFFUN(exprs(X), clusters, pval.type="all")
    REFFUN(exprs(X), clusters, direction="up", pval.type="all")
    REFFUN(exprs(X), clusters, direction="down", pval.type="all")

    # Checking what happens if one of the groups has only one element.
    re.clust <- clust$cluster
    re.clust[1] <- 4
    re.clust <- factor(re.clust)
    REFFUN(exprs(X), re.clust)
    REFFUN(exprs(X), re.clust, pval.type="all")

    # Checking what happens if two of the groups have only one element.
    re.clust <- clust$cluster
    re.clust[1:2] <- 4:5
    re.clust <- factor(re.clust)
    REFFUN(exprs(X), re.clust)
    REFFUN(exprs(X), re.clust, pval.type="all")
    
    # Checking what happens if there is an empty level.
    re.clusters <- clusters
    levels(re.clusters) <- 1:4

    out <- findMarkers(exprs(X), re.clusters)
    ref <- findMarkers(exprs(X), clusters)
    for (g in names(ref)) {
        current <- ref[[g]]
        counter <- out[[g]]
        expect_equal(current, counter[,colnames(current)])
        expect_true(all(is.na(counter[,"logFC.4"])))
        expect_true(all(is.na(out[["4"]][,paste0("logFC.", g)])))
    }
})

test_that("findmarkers works with a log-fold change threshold", {
    clusters <- rep(1:2, ncells/2)             
    exprs(X) <- exprs(X) + outer(rnorm(ngenes, sd=1), as.integer(clusters))

    # Log-fold change with direction="any".
    out <- findMarkers(exprs(X), clusters, lfc=1, pval.type="all")
    N <- ncol(out[[1]])
    for (xx in out) {
        tab <- table(xx$IUT.p < 0.5, abs(xx[[N]]) > 1)
        expect_identical(sum(diag(tab)), nrow(xx))
    }
    out <- findMarkers(exprs(X), clusters, lfc=50, pval.type="all")
    for (xx in out) {
        expect_equal(xx$IUT.p, rep(1, ngenes)) 
    }

    # Log-fold change with direction="up".
    out <- findMarkers(exprs(X), clusters, lfc=1, direction="up", pval.type="all")
    for (xx in out) {
        tab <- table(xx$IUT.p < 0.5, xx[[N]] > 1)
        expect_identical(sum(diag(tab)), nrow(xx))
    }
    out <- findMarkers(exprs(X), clusters, lfc=50, direction="up", pval.type="all")
    for (xx in out) {
        expect_equal(xx$IUT.p, rep(1, ngenes)) 
    }

    # Log-fold change with direction="down".
    out <- findMarkers(exprs(X), clusters, lfc=1, direction="down", pval.type="all")
    for (xx in out) {
        tab <- table(xx$IUT.p < 0.5, xx[[N]] < -1)
        expect_identical(sum(diag(tab)), nrow(xx))
    }
    out <- findMarkers(exprs(X), clusters, lfc=50, direction="down", pval.type="all")
    for (xx in out) {
        expect_equal(xx$IUT.p, rep(1, ngenes)) 
    }
})

###############################
# Checking that the blocking runs. Unfortunately, the output of unblocked
# findMarkers cannot be easily used to evaluated the correctness of the blocked
# version, as the results are not simple sums. Instead, we have to set up 
# scenarios where the output is easier to check.

set.seed(70000011)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

clust <- kmeans(t(exprs(X)), centers=3)
clusters <- as.factor(clust$cluster)
blocked <- sample(3, ncells, replace=TRUE)

test_that("findMarkers runs properly with blocking: basic checks", {
    # All values should be the same if one block contains only one group.
    re.clusters <- clusters
    re.clusters[blocked==1] <- 1
    output <- findMarkers(X, re.clusters, blocked)
    keep <- blocked!=1
    ref <- findMarkers(X[,keep], re.clusters[keep], blocked[keep])
    expect_equal(ref, output)

    # If one block contains only one sample from each group, gene ordering and p-values 
    # should be the same as an unblocked analysis with the other block.
    re.clusters <- clusters
    re.blocked <- rep(1, length(clusters))
    re.clusters[1:3] <- 1:3
    re.blocked[1:3] <- 2

    for (type in c("any", "all")) {
        if (type=="any") {
            direction <- c("any", "up", "down")
        } else {
            direction <- "any"
        }

        for (d in direction) { 
            output <- findMarkers(X, re.clusters, block=re.blocked, direction=d, pval.type=type)
            ref <- findMarkers(X[,-(1:3)], re.clusters[-(1:3)], direction=d, pval.type=type)
            expect_equal(names(output), names(ref))

            for (i in names(output)) { 
                expect_equal(rownames(output[[i]]), rownames(ref[[i]]))
                expect_equal(output[[i]]$FDR, ref[[i]]$FDR)
                expect_equal(output[[i]]$Top, ref[[i]]$Top)
                expect_equal(output[[i]]$IUT.p, ref[[i]]$IUT.p)
            }
        }
    }

    # Ordering and log-fold changes should be the same if the block is duplicated.
    clusters2 <- rep(clusters, 2)
    blocked2 <- rep(1:2, each=length(clusters))
    X2 <- BiocGenerics::cbind(X, X)

    for (type in c("any", "all")) {
        if (type=="any") {
            direction <- c("any", "up", "down")
        } else {
            direction <- "any"
        }

        for (d in direction) { 
            output <- findMarkers(X2, clusters2, block=blocked2, direction=d, pval.type=type)
            ref <- findMarkers(X, clusters, direction=d, pval.type=type)
            expect_equal(names(output), names(ref))

            for (i in names(output)) { 
                expect_equal(rownames(output[[i]]), rownames(ref[[i]]))
                expect_equal(output[[i]]$Top, ref[[i]]$Top)
                expect_equal(output[[i]][,-(1:2)], ref[[i]][,-(1:2)])
            }
        }
    }
})

set.seed(70000012)
test_that("findMarkers runs properly with blocking: combining checks", {
    # Checking that the combined log-fold changes follow the expected ratios,
    # for a multi-cluster comparison with three blocking levels but only 
    # two unique "X". Note we change the logFCs without altering the variance.
    clusters3 <- rep(clusters, 3)
    blocked3 <- rep(1:3, each=length(clusters))
    new.X <- X
    exprs(new.X) <- exprs(new.X) + outer(runif(ngenes), as.integer(clusters)) 
    X3 <- BiocGenerics::cbind(X, X, new.X)

    output <- findMarkers(X3, clusters3, block=blocked3)
    ref1 <- findMarkers(X, clusters)
    ref2 <- findMarkers(new.X, clusters)

    for (i in names(output)) {
        cur.out <- output[[i]][rownames(X),-(1:2)]
        r.out1 <- ref1[[i]][rownames(X),-(1:2)]
        r.out2 <- ref2[[i]][rownames(X),-(1:2)]
        for (j in colnames(cur.out)) { 
            expect_equal(cur.out[,j], (r.out1[,j]*2+r.out2[,j])/3)
        }
    }

    # Checking that the combined p-values follow the expected ratios
    # (only setting 1 pair of groups with pval.type="all", so that 
    # "IUT.p" can be directly extracted as the p-value). Note that 
    # the log-fold change shift has to be small to avoid breaking qnorm.
    new.clusters <- sample(2, ncol(X), replace=TRUE)
    new.clusters3 <- rep(new.clusters, 3)
    new.X <- X
    exprs(new.X) <- exprs(new.X) + outer(rnorm(ngenes, sd=0.01), as.integer(new.clusters))
    new.X3 <- BiocGenerics::cbind(X, X, new.X)

    for (d in c("up", "down")) { 
        output <- findMarkers(new.X3, new.clusters3, block=blocked3, pval.type="all", direction=d)
        ref1 <- findMarkers(X, new.clusters, pval.type="all", direction=d)
        ref2 <- findMarkers(new.X, new.clusters, pval.type="all", direction=d)

        for (i in names(output)) {
            cur.out <- output[[i]][rownames(X),"IUT.p"]
            r.out1 <- ref1[[i]][rownames(X),"IUT.p"]
            r.out2 <- ref2[[i]][rownames(X),"IUT.p"]
            expect_equal(qnorm(cur.out), (qnorm(r.out1)*2+qnorm(r.out2))/3)
        }
    }
})

set.seed(70000012)
test_that("findMarkers runs properly with blocking: missing/singleton group checks", {
    # Log-fold changes behave sensibly when one block is missing a group.
    re.clusters <- clusters
    re.blocked <- blocked
    re.clusters[re.blocked==1 & re.clusters==1] <- 2

    output <- findMarkers(X, re.clusters, block=re.blocked)
    keep <- re.blocked!=1
    output2 <- findMarkers(X[,keep], re.clusters[keep], block=re.blocked[keep])
    expect_equal(output[["1"]], output2[["1"]])

    # Log-fold changes are just block averages if there are no residual d.f. anywhere.
    subX <- X[,1:9]
    subclust <- rep(1:3, each=3)
    subblock <- rep(1:3, 3)

    output <- findMarkers(subX, subclust, subblock)
    for (i in names(output)) {
        expect_true(all(is.na(output[[i]]$Top)))
        for (j in setdiff(names(output), i)) {
            cur.lfc <- output[[i]][,paste0("logFC.", j)]
            expect_equal(cur.lfc, unname(rowMeans(exprs(subX)[,subclust==i]) - rowMeans(exprs(subX)[,subclust==j])))
        }
    }

    # Manufacturing a situation where the variance is easily calculated, to test the 
    # inference of the standard error when there is no d.f. in one block.
    stuff <- matrix(rnorm(ngenes*ncells), nrow=ngenes)
    rownames(stuff) <- 1:ngenes
    re.clusters <- rep(1:2, length.out=ncells)
    re.block <- rep(1, ncells)
    re.block[1:2] <- 2

    output <- findMarkers(stuff, re.clusters, re.block) 
    in.11 <- stuff[,re.block==1 & re.clusters==1]
    in.12 <- stuff[,re.block==1 & re.clusters==2]
    w.b2 <- 1/2
    w.b1 <- 1/(1/ncol(in.11) + 1/ncol(in.12))
    expect_equal(output[[2]][rownames(stuff),]$logFC.1, 
                 unname(((stuff[,2]-stuff[,1]) * w.b2 + 
                         (rowMeans(in.12) - rowMeans(in.11)) * w.b1)/(w.b2 + w.b1))
                 )
    
    re.clusters[] <- 2
    re.clusters[1:2] <- 1:2
    output <- findMarkers(stuff, re.clusters, re.block) 
    expect_equal(output[[2]]$logFC.1, unname(stuff[,2]-stuff[,1]))
    
    # We correctly get NA values if block is confounded with group.
    # Importantly, these NA values should not affect anything else.
    stuff <- matrix(rnorm(ngenes*ncells), nrow=ngenes)
    rownames(stuff) <- 1:ngenes
    re.clusters <- as.character(clusters)
    re.block <- blocked
    re.clusters[re.block==1] <- "A"

    output <- findMarkers(stuff, re.clusters, re.block) 
    keep <- re.block != 1
    ref <- findMarkers(stuff[,keep], re.clusters[keep], re.block[keep]) 

    for (x in 1:3) {
        cur.out <- output[[x]]
        expect_true(all(is.na(cur.out$logFC.A)))
        expect_equal(cur.out[,setdiff(colnames(cur.out), "logFC.A")], ref[[x]])
        expect_true(all(is.na(output[["A"]][,paste0("logFC.", x)])))
    } 
})

###############################

# Setting up a reference function for the log-fold changes ONLY.
# This is because computing the p-values would require me to just
# copy the R code over, which defeats the purpose.
library(limma)
LMFUN <- function(y, design, clust.vals, output) { 
    lfit <- lmFit(y, design)
    for (host in clust.vals) {
        collected.lfc <- collected.p <- list()
        
        for (target in setdiff(clust.vals, host)) {
            con <- setNames(numeric(ncol(design)), colnames(design))
            con[[host]] <- 1
            con[[target]] <- -1

            fit2 <- contrasts.fit(lfit, con)
            fit2 <- eBayes(fit2)
            res <- topTable(fit2, n=Inf, sort.by="none")
            collected.lfc[[paste0("logFC.", target)]] <- res$logFC
        }
        
        combined.lfc <- do.call(cbind, collected.lfc)
        rownames(combined.lfc) <- rownames(y)
        cur.out <- output[[host]]
        expect_equal(as.matrix(cur.out[rownames(y),-(1:2)]), combined.lfc)
    }  
    return(TRUE)
}

set.seed(7000002)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

test_that("findMarkers works properly with a design matrix", {
    clust <- kmeans(t(exprs(X)), centers=3)
    clusters <- as.factor(clust$cluster)

    # Trying with an additive blocking factor.
    block <- factor(sample(2, ncol(X), replace=TRUE))
    design <- model.matrix(~0 + clusters + block)
    colnames(design) <- c(levels(clusters), "block2")

    out.des <- findMarkers(X, clusters=clust$cluster, design=model.matrix(~block))
    LMFUN(exprs(X), design, levels(clusters), out.des)

    # Trying with a real-value blocking factor.
    covariates <- runif(ncells)
    design <- model.matrix(~0 + clusters + covariates)
    colnames(design) <- c(levels(clusters), "covariate")

    out.des <- findMarkers(X, clusters=clust$cluster, design=model.matrix(~covariates))
    LMFUN(exprs(X), design, levels(clusters), out.des)

    # Checking for correct error upon having no residual d.f. in the design matrix.
    design0 <- matrix(1, ncol(X), 1)
    expect_error(findMarkers(X, clusters=1:ncol(X), design=design0), "no residual d.f.")

    # Checking that block= takes priority.
    outdes <- findMarkers(X, clusters=clust$cluster, block=block, design=model.matrix(~block))
    outref <- findMarkers(X, clusters=clust$cluster, block=block)
    expect_equal(outdes, outref)
})

test_that("findMarkers with a design matrix responds to log-fold change thresholds", {
    clusters <- rep(1:2, ncells/2)             
    exprs(X) <- exprs(X) + outer(rnorm(ngenes, sd=1), as.integer(clusters))
    block <- factor(sample(2, ncol(X), replace=TRUE))
    design <- model.matrix(~block)

    # Log-fold change with direction="any".
    out <- findMarkers(exprs(X), clusters, lfc=1, pval.type="all", design=design)
    N <- ncol(out[[1]])
    for (xx in out) {
        tab <- table(xx$IUT.p < 0.5, abs(xx[[N]]) > 1)
        expect_identical(sum(diag(tab)), nrow(xx))
    }
    out <- findMarkers(exprs(X), clusters, lfc=50, pval.type="all", design=design)
    for (xx in out) {
        expect_equal(xx$IUT.p, rep(1, ngenes)) 
    }

    # Log-fold change with direction="up".
    out <- findMarkers(exprs(X), clusters, lfc=1, direction="up", pval.type="all", design=design)
    for (xx in out) {
        tab <- table(xx$IUT.p < 0.5, xx[[N]] > 1)
        expect_identical(sum(diag(tab)), nrow(xx))
    }
    out <- findMarkers(exprs(X), clusters, lfc=50, direction="up", pval.type="all", design=design)
    for (xx in out) {
        expect_equal(xx$IUT.p, rep(1, ngenes)) 
    }

    # Log-fold change with direction="down".
    out <- findMarkers(exprs(X), clusters, lfc=1, direction="down", pval.type="all", design=design)
    for (xx in out) {
        tab <- table(xx$IUT.p < 0.5, xx[[N]] < -1)
        expect_identical(sum(diag(tab)), nrow(xx))
    }
    out <- findMarkers(exprs(X), clusters, lfc=50, direction="down", pval.type="all", design=design)
    for (xx in out) {
        expect_equal(xx$IUT.p, rep(1, ngenes)) 
    }
})

###############################

set.seed(7000003)
test_that("findMarkers works correctly with subsetting and spikes", {   
    # Works with an SCE object.
    clust <- kmeans(t(exprs(X)), centers=3)
    out <- findMarkers(X, clusters=clust$cluster)
    out2 <- findMarkers(exprs(X), clusters=clust$cluster)
    expect_identical(out, out2)

    # Works with subsetting and spikes.
    out <- findMarkers(X, clusters=clust$cluster, subset.row=100:1)
    out2 <- findMarkers(X[100:1,], clusters=clust$cluster)
    expect_identical(out, out2)
    
    isSpike(X, "ERCC") <- 1:100
    out <- findMarkers(X, clusters=clust$cluster)
    out2 <- findMarkers(exprs(X)[-(1:100),], clusters=clust$cluster)
    expect_identical(out, out2)

    # Repeating with a design matrix, to check that subsetting works in both branches for coefficient calculation.
    block <- factor(sample(2, ncol(X), replace=TRUE))
    out.des <- findMarkers(exprs(X), clusters=clust$cluster, design=model.matrix(~block), subset.row=100:1)
    out.des2 <- findMarkers(exprs(X)[100:1,,drop=FALSE], clusters=clust$cluster, design=model.matrix(~block))
    expect_identical(out.des, out.des2)
})

test_that("log-transformed p-values in findMarkers works correctly", {
    # Checking that the log-transformed argument works.
    default <- findMarkers(exprs(X), clusters)
    logged <- findMarkers(exprs(X), clusters, log.p=TRUE)
    expect_identical(names(default), names(logged))
    for (i in names(default)) {
        expect_identical(rownames(default[[i]]), rownames(logged[[i]]))
        expect_equal(log(default[[i]]$FDR), logged[[i]]$log.FDR)
    }

    default <- findMarkers(exprs(X), clusters, direction="up")
    logged <- findMarkers(exprs(X), clusters, log.p=TRUE, direction="up")
    expect_identical(names(default), names(logged))
    for (i in names(default)) {
        expect_identical(rownames(default[[i]]), rownames(logged[[i]]))
        expect_equal(log(default[[i]]$FDR), logged[[i]]$log.FDR)
    }

    default <- findMarkers(exprs(X), clusters, pval.type="all")
    logged <- findMarkers(exprs(X), clusters, log.p=TRUE, pval.type="all")
    expect_identical(names(default), names(logged))
    for (i in names(default)) {
        expect_identical(rownames(default[[i]]), rownames(logged[[i]]))
        expect_equal(log(default[[i]]$IUT.p), logged[[i]]$log.IUT.p)
        expect_equal(log(default[[i]]$FDR), logged[[i]]$log.FDR)
    }

    # Checking that the log-transforming works for full.stats=TRUE.
    default <- findMarkers(exprs(X), clusters, full.stats=TRUE)
    logged <- findMarkers(exprs(X), clusters, log.p=TRUE, full.stats=TRUE)
    expect_identical(names(default), names(logged))
    for (i in names(default)) {
        curdefault <- default[[i]]
        curlogged <- logged[[i]]
        expect_identical(rownames(curdefault), rownames(curlogged))
        for (j in grep("^stats", colnames(curdefault))) {
            expect_equal(log(default[[i]][[j]]$FDR), logged[[i]][[j]]$log.FDR)
            expect_equal(log(default[[i]][[j]]$p.value), logged[[i]][[j]]$log.p.value)
        }
    }
    
    # Also for linear models.
    block <- factor(sample(2, ncol(X), replace=TRUE))
    design <- model.matrix(~block)
    default <- findMarkers(exprs(X), clusters, design=design)
    logged <- findMarkers(exprs(X), clusters, log.p=TRUE, design=design)
    expect_identical(names(default), names(logged))
    for (i in names(default)) {
        expect_identical(rownames(default[[i]]), rownames(logged[[i]]))
        expect_equal(log(default[[i]]$FDR), logged[[i]]$log.FDR)
    }
})

# Checking consistency upon silly inputs.

test_that("findMarkers behaves sensibly with silly inputs", {
    # No genes.
    out <- findMarkers(exprs(X)[0,], clusters=rep(1:3, length.out=ncol(X)))
    for (i in out) { expect_identical(nrow(i), 0L) }    
    
    # No cells.
    out <- findMarkers(exprs(X)[,0], clusters=integer(0))
    expect_identical(length(out), 0L)    
    expect_error(findMarkers(exprs(X)[,0], clusters=integer(0), design=matrix(0,0,0)), 
                 "contrasts can be applied only to factors with 2 or more levels")

    # Mismatch in dimensions.
    expect_error(findMarkers(exprs(X), clusters=1), "length of 'clusters' does not equal 'ncol(x)'", fixed=TRUE)
    dummy <- sample(2, ncol(X), replace=TRUE)
    expect_error(findMarkers(exprs(X), clusters=dummy, block=1), "length of 'block' does not equal 'ncol(x)'", fixed=TRUE)
    expect_error(findMarkers(exprs(X), clusters=dummy, design=matrix(0,0,0)), "'nrow(design)' is not equal to 'ncol(x)'", fixed=TRUE)

    # Checking that we get non-NA p-values with no variance.
    clusters <- rep(1:2, each=ncells/2)
    stuff <- matrix(clusters, ngenes, ncells, byrow=TRUE)
    out <- findMarkers(stuff, clusters)
    expect_true(all(out[[1]]$FDR < 0.01))
    expect_true(all(out[[2]]$FDR < 0.01))
    expect_equal(out[[1]]$logFC.2, rep(-1, ngenes))
    expect_equal(out[[2]]$logFC.1, rep(1, ngenes))
})

