#########################################
# Loading packages
#########################################

library(MERINGUE)
library(crawdad)
library(tidyverse)
library(ggplot2)

#########################################
# Pre-processing
#########################################

loadingData <- function(filename1, filename2,
                        xmin = 5000, xmax = 6000,
                        ymin = 5000, ymax = 6000) {
  
  counts <- fread(filename1)
  meta   <- fread(filename2, select = c('center_x','center_y'))
  
  colnames(meta) <- c('x', 'y')
  
  ## Plot 1: The entire feature plot
  meta %>%
    ggplot() +
    geom_point(aes(x, y), colour = 'lightgrey', size = .0001)
  
  meta <- meta %>%
    filter(x > xmin, x < xmax, y > ymin, y < ymax)
  
  ## Plot 2: The cropped feature plot
  meta %>%
    ggplot() +
    geom_point(aes(x, y), colour = 'lightpink', size = .5)
  
  counts<- counts %>%
    filter(rownames(counts) %in% rownames(meta))

  return(list(counts, meta))
}

preProcessing <- function(cd, pos, 
                          reads    = 200,
                          lib.size = 200) {
  counts <- cleanCounts(counts       = t(cd), 
                        min.reads    = reads, 
                        min.lib.size = lib.size, 
                        plot         = TRUE,
                        verbose      = TRUE)
  
  pos <- pos %>%
    filter(rownames(pos) %in% col(counts))
  
  ## Plot 3: The cropped + filtered feature plot
  pos %>% ggplot() +
    geom_point(aes(x,y), colour = 'lightpink', size = .5)
  
  # Normalize count
  mat <- normalizeCounts(counts = counts, 
                         normFactor = 1,
                         log=FALSE,
                         verbose=TRUE)
  
  return(list(pos, mat))
}

lung.data <- loadingData("HumanLungCancerPatient1_cell_by_gene.csv",
                         "HumanLungCancerPatient1_cell_metadata.csv")
cd  <- lung.data[[1]]
pos <- lung.data[[2]]
process <- preProcessing(cd, pos)
pos <- process[[1]]
mat <- process[[2]]

#########################################
# Clustering
#########################################

# PCA Dimension reduction on log10 CPM expression values
pcs.info <- prcomp(t(log10(as.matrix(mat)+1)), center=TRUE)

pcs.sdev <- pcs.info$sdev

ggplot() + geom_point(aes(1:20,pcs.sdev[1:20]))

pcs.sdev[1:20]

# manually decided
nPcs <- 15

pcs <- pcs.info$x[,1:nPcs]

# 2D embedding by tSNE for visualization
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y

# both vars have the same row name
rownames(emb) <- rownames(pcs) <- 1:dim(pcs)[1]

# Graph-based cluster detection
k <- 30 # linking 30 closest cells to cell of interest
com <- getClusters(pcs, k, weight=TRUE)

# Manually annotate identified clusters with cell-types
annot <- as.character(com)
names(annot) <- names(com)
annot <- as.factor(annot)

# Plot
par(mfrow=c(1,2), mar=rep(1,4))
plotEmbedding(emb, groups=annot, 
              show.legend=TRUE, xlab=NA, ylab=NA,
              verbose=FALSE)
plotEmbedding(pos, groups=annot, 
              cex=1, xlab=NA, ylab=NA,
              verbose=FALSE)

#########################################
# SVG Identification
#########################################

colnames(mat) <- 1:dim(mat)[2]
dg <- getDifferentialGenes(as.matrix(mat), annot)


dg.sig <- lapply(dg, function(x) {
  x <- x[x$p.adj < 0.05,]
  x <- na.omit(x)
  x <- x[x$highest,]
  rownames(x)
})
print(lapply(dg.sig, length))

dg.genes <- unlist(dg.sig)
ggroup <- unlist(lapply(1:length(dg.sig), function(i) { 
  rep(names(dg.sig)[i], length(dg.sig[[i]]))
}))
names(ggroup) <- dg.genes
ggroup <- factor(ggroup)

# Plot
ccol <- rainbow(length(levels(annot)))[annot]
names(ccol) <- names(annot) # column colors
gcol <- rainbow(length(levels(ggroup)), v=0.5)[ggroup]
names(gcol) <- names(ggroup) # row colors

m <- as.matrix(mat[dg.genes, names(sort(annot))])
m <- winsorize(t(scale(t(m))))

heatmap(m, scale="none", 
        Colv=NA, Rowv=NA, labRow=NA, labCol=NA,
        ColSideColors=ccol[colnames(m)],
        RowSideColors=gcol[rownames(m)],
        col=colorRampPalette(c('blue', 'white', 'red'))(100)
)


# spatially informed analysis

# Get neighbor-relationships
w <- getSpatialNeighbors(pos, verbose = TRUE, filterDist = 4)
# 0 1
# 63489077 47764

plotNetwork(pos, w)

# Identify significantly spatially auto-correlated genes
start_time <- Sys.time()
I <- getSpatialPatterns(mat, w)
end_time <- Sys.time()
print(end_time - start_time) 
# Time difference of 20.99901  mins

results.filter <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)

# Number of significantly auto-correlated genes: 10
# Compute spatial cross correlation matrix
# ...driven by > 397.55 cells: 10

scc <- spatialCrossCorMatrix(mat = as.matrix(mat[results.filter,]), 
                             weight = w)
# Identify primary patterns
par(mfrow=c(2,2), mar=rep(2,4))
ggroup <- groupSigSpatialPatterns(pos = pos, 
                                  mat = as.matrix(mat[results.filter,]), 
                                  scc = scc, 
                                  power = 1, 
                                  hclustMethod = 'ward.D', 
                                  deepSplit = 2,
                                  zlim=c(-1.5,1.5))
# groups
# 1 2 3 4
# 4 2 2 2

#########################################
# cell-type colocalization
#########################################

df <- data.frame(x = pos$x, y = pos$y, celltypes = annot)

df %>% ggplot() + geom_point(aes(x, y, color= as.numeric(annot)))

## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = df[,c("x","y")], celltypes = df$celltypes)

## define the scales to analyze the data
scales <- c(100, 200, 300, 400, 500, 600, 700)

shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = 14,
                                            seed = 1,
                                            verbose = TRUE)

saveRDS(shuffle.list, "myshufflelist.RDS")

## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells,
                               dist = 4,
                               shuffle.list = shuffle.list,
                               ncores = 14,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)

zsig <- correctZBonferroni(dat)

## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
