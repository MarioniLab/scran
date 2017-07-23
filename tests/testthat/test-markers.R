# This tests the findMarkers function.
# require(scran); require(testthat); source("test-markers.R")

set.seed(70000)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

# Setting up a reference function.
library(limma)
REFFUN <- function(y, design, clust.vals, output, pval.type="any", direction="any", min.mean=0.1) { 
    lfit <- lmFit(y, design)
    all.means <- rowMeans(y)
    do.solo <- is.null(min.mean) || all(all.means >= min.mean) || !any(all.means >= min.mean)

    for (host in clust.vals) {
        collected.lfc <- collected.p <- list()
        
        for (target in setdiff(clust.vals, host)) {
            con <- setNames(numeric(ncol(design)), colnames(design))
            con[[host]] <- 1
            con[[target]] <- -1
            
            fit2 <- contrasts.fit(lfit, con)
            if (do.solo) {
                fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
                res <- topTable(fit2, n=Inf, sort.by="none")
            } else {
                higher <- all.means >= min.mean
                fit2a <- fit2[higher,]
                fit2b <- fit2[!higher,]
                fit2a <- eBayes(fit2a, trend=TRUE, robust=TRUE)
                fit2b <- eBayes(fit2b, trend=TRUE, robust=TRUE)

                resa <- topTable(fit2a, n=Inf, sort.by="none")
                resb <- topTable(fit2b, n=Inf, sort.by="none")
                res <- rbind(resa, resb)
                res <- res[order(c(which(higher), which(!higher))),]                
            }

            if (direction=="up") {
                is.up <- res$logFC > 0
                res$P.Value[is.up] <- res$P.Value[is.up]/2
                res$P.Value[!is.up] <- 1-res$P.Value[!is.up]/2
            } else if (direction=="down") {
                is.down <- res$logFC < 0
                res$P.Value[is.down] <- res$P.Value[is.down]/2
                res$P.Value[!is.down] <- 1-res$P.Value[!is.down]/2
            }
            
            collected.lfc[[paste0("logFC.", target)]] <- res$logFC
            collected.p[[target]] <- res$P.Value
        }
        
        # Checking the value of 'top'.
        all.ranks <- lapply(collected.p, rank, ties.method="first")
        cur.out <- output[[host]]
        for (i in seq(1, ngenes, length.out=51)) {
            current <- rownames(X)[unique(unlist(lapply(all.ranks, FUN=function(x) { which(x <= i) })))]
            expect_identical(sort(current), sort(cur.out$Gene[cur.out$Top <= i]))
        }
       
        reordered <- order(do.call(pmin, all.ranks), do.call(pmin, collected.p))
        expect_identical(cur.out$Gene, rownames(y)[reordered])
        
        # Checking the log-fold changes.
        for (target in names(collected.lfc)) {
            expect_equal(collected.lfc[[target]][reordered], cur.out[[target]])
        }
        
        # Checking the FDR, after Simes' or with the IUT.
        combined.p <- do.call(cbind, collected.p)
        if (pval.type=="any") { 
            pval <- apply(combined.p, 1, FUN=function(x) { min(p.adjust(x, method="BH")) })
        } else {
            pval <- apply(combined.p, 1, FUN=max)
        }
        adj.p <- p.adjust(pval, method="BH")
        expect_equal(adj.p[reordered], cur.out$FDR)
    }  
    return(TRUE)
}

clust <- kmeans(t(exprs(X)), centers=3)
out <- findMarkers(X, clusters=clust$cluster)

clusters <- as.factor(clust$cluster)
design <- model.matrix(~0 + clusters)
colnames(design) <- levels(clusters)
REFFUN(exprs(X), design, levels(clusters), out)

# Checking that the IUT calculation is correct.

out <- findMarkers(X, clusters=clust$cluster, pval.type="all")
REFFUN(exprs(X), design, levels(clusters), out, pval.type="all")

# Checking that the directional calculations are correct.

out <- findMarkers(X, clusters=clust$cluster, direction="up")
REFFUN(exprs(X), design, levels(clusters), out, direction="up")

out <- findMarkers(X, clusters=clust$cluster, direction="down")
REFFUN(exprs(X), design, levels(clusters), out, direction="down")

# Checking how it behaves with a design matrix.

clust <- kmeans(t(exprs(X)), centers=3)
block <- factor(sample(2, ncol(X), replace=TRUE))
out.des <- findMarkers(X, clusters=clust$cluster, design=model.matrix(~block))

clusters <- as.factor(clust$cluster)
design <- model.matrix(~0 + clusters + block)
colnames(design) <- c(levels(clusters), "block2")
REFFUN(exprs(X), design, levels(clusters), out.des)

# Checking that the results are the same regardless of how I provide inputs.

out <- findMarkers(X, clusters=clust$cluster)
out2 <- findMarkers(exprs(X), clusters=clust$cluster)
expect_identical(out, out2)

out <- findMarkers(X, clusters=clust$cluster, subset.row=100:1)
out2 <- findMarkers(X[100:1,], clusters=clust$cluster)
expect_identical(out, out2)

X <- calculateQCMetrics(X, feature_controls=list(ERCC=1:100))
setSpike(X) <- "ERCC"
out <- findMarkers(X, clusters=clust$cluster)
out2 <- findMarkers(exprs(X)[-(1:100),], clusters=clust$cluster)
expect_identical(out, out2)

# Repeating with non-infinite d.f. to check shrinkage.

set.seed(700001)
test_that("findMarkers works with non-infinite prior d.f.", {
    s2 <- 10/rchisq(1000, df=10)
    y <- matrix(rnorm(1000*200, sd=sqrt(s2)), nrow=1000)
    rownames(y) <- paste0("X", seq_len(1000))
    clusters <- factor(sample(4, 200, replace=TRUE))
    out <- findMarkers(y, clusters=clusters)
    
    design <- model.matrix(~0 + clusters)
    colnames(design) <- levels(clusters)
    REFFUN(y, design, levels(clusters), out)

    # Trying again with fewer prior d.f.    
    s2 <- 5/rchisq(1000, df=5)
    y <- matrix(rnorm(1000*100, sd=sqrt(s2)), nrow=1000)
    rownames(y) <- paste0("X", seq_len(1000))
    clusters <- factor(sample(5, 100, replace=TRUE))
    out <- findMarkers(y, clusters=clusters)
    
    design <- model.matrix(~0 + clusters)
    colnames(design) <- levels(clusters)
    REFFUN(y, design, levels(clusters), out)
})

# Testing the min.mean setting.

set.seed(700002)
test_that("findMarkers works with a variety of minimum means", {
    s2 <- 10/rchisq(1000, df=10)
    y <- matrix(rnorm(1000*200, sd=sqrt(s2)), nrow=1000)
    rownames(y) <- paste0("X", seq_len(1000))
    clusters <- factor(sample(3, 200, replace=TRUE))
    
    design <- model.matrix(~0 + clusters)
    colnames(design) <- levels(clusters)
    out <- findMarkers(y, clusters=clusters, min.mean=0)
    REFFUN(y, design, levels(clusters), out, min.mean=0)

    out <- findMarkers(y, clusters=clusters, min.mean=NULL)
    REFFUN(y, design, levels(clusters), out, min.mean=NULL)
})

# Checking consistency upon silly inputs.

expect_error(findMarkers(exprs(X)[0,], clusters=clust$cluster), "var is empty")
expect_error(findMarkers(exprs(X)[,0], clusters=integer(0)), "contrasts can be applied only to factors with 2 or more levels")

