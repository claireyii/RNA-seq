load("../data/se_paul.rda")
library(SingleCellExperiment)
library(scater)

sim <- SingleCellExperiment(se)
counts <- counts(assay(sim))
sim <- SingleCellExperiment(assays = List(counts = counts))

col_dat <- colData(sim)

sim <- logNormCounts(sim)
sim <- runPCA(sim)
sim <- runUMAP(sim)

# How to construct clusters for Slingshot?
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sim)$GMM <- cl1 # sim did not have cell types in colData
plotReducedDim(sim, dimred = "PCA", colour_by = "GMM") 

# Slingshot
library(slingshot)
library(tradeSeq)

cl1 <- Mclust(rd1)$classification
colData(sim)$GMM <- cl1
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')

plot(reducedDims(sim)$PCA, col = brewer.pal(9,'Set1')[sim$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')

lin1 <- getLineages(rd1, cl, start.clus= '1', end.clus = '2')
plot(rd1, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(lin1, lwd = 3, col = 'black', show.constraints = TRUE)