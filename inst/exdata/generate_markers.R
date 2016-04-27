library(scran)
dir.create("temp")

# MOUSE (relies on semi-public data):

zipfile <- "temp/current.zip"
download.file("https://github.com/PMBio/cyclone/archive/c982e7388d8e49e1459055504313f87bb3eb0ceb.zip", zipfile)
out <- unzip(zipfile, exdir="temp")

pdata <- out[grep("pairs_functions.RData$", out)]
load(pdata)

all.pairs <- sandbag(is.G1=id.G1, is.S=id.S, is.G2M=id.G2M, training.data[genes.training,], fraction=0.5)
saveRDS(file="mouse_cycle_markers.rds", all.pairs)

rm(list=ls())

# HUMAN:

library(GEOquery)
out <- getGEOSuppFiles("GSE64016", baseDir="temp", makeDirectory=FALSE)
count.file <- "temp/GSE64016_H1andFUCCI_normalized_EC.csv.gz"
hs.counts <- read.csv(count.file, header=TRUE, row.names=1)
hs.G1 <- grepl("G1", colnames(hs.counts))
hs.S <- grepl("S", colnames(hs.counts))
hs.G2 <- grepl("G2", colnames(hs.counts))

library(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, keytype="SYMBOL", key=rownames(hs.counts), column="ENSEMBL")
anno <- anno[!is.na(anno$ENSEMBL),]
m <- match(anno$SYMBOL, rownames(hs.counts))
hs.counts2 <- as.matrix(hs.counts[m,])
rownames(hs.counts2) <- anno$ENSEMBL
hs.cycle <- select(org.Hs.eg.db, keytype="GOALL", key="GO:0007049", column="ENSEMBL")
hs.training <- rownames(hs.counts2) %in% hs.cycle$ENSEMBL

all.pairs <- sandbag(is.G1=hs.G1, is.S=hs.S, is.G2M=hs.G2, hs.counts2[hs.training,], fraction=0.5)
saveRDS(file="human_cycle_markers.rds", all.pairs)

# Cleaning up:

unlink("temp", recursive=TRUE)
