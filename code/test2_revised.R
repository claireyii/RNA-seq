load("../data/se_paul.rda")
library(SingleCellExperiment)
library(scater)

# sim <- SingleCellExperiment(se)
# counts <- counts(assay(sim))
# sim <- SingleCellExperiment(assays = List(counts = counts))

sim <- se

counts(sim) <- as.matrix(counts(sim))

col_dat <- colData(sim)

sim <- logNormCounts(sim)
sim <- runPCA(sim)
sim <- runUMAP(sim)

# How to construct clusters for Slingshot?
rd1 <- reducedDim(sim, "PCA")
library(mclust, quietly = TRUE)
cl1 <- kmeans(rd1, 7)$cluster
sim$kmeans <- cl1 # sim did not have cell types in colData
plotReducedDim(sim, dimred = "UMAP", colour_by = "kmeans") 
plotReducedDim(sim, dimred = "UMAP", colour_by = "cell_type2") 

# Slingshot
library(slingshot)
library(tradeSeq)

# cl1 <- Mclust(rd1)$classification
# colData(sim)$GMM <- cl1
library(RColorBrewer)
plot(rd1, col = brewer.pal(9, "Set1")[cl1], pch = 16, asp = 1)

sim <- slingshot(sim, clusterLabels = "kmeans", reducedDim = "UMAP",
                 start.clus = 3)

plot(reducedDims(sim)$UMAP, col = brewer.pal(9, "Set1")[sim$kmeans], pch = 16, asp = 1)
lines(SlingshotDataSet(sim), lwd = 2, col = "black")

lin1 <- getLineages(rd1, cl1, start.clus = "1", end.clus = "2")
plot(rd1, col = brewer.pal(9, "Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3, col = "black", show.constraints = TRUE)

# tradeSeq 
crv <- getCurves(lin1)
sce <- fitGAM(counts = counts(sim), sds = crv)
saveRDS(sce, file = "fitGAM_output.rds")
prediction <- predictSmooth(sce, gene = 1:10, nPoints = 40)


# A deeper look at the curves
## Extract the curve objects
sds <- SlingshotDataSet(sim)
## Look at the principal curve part
?principal_curve
## espeically this object
sds@curves$curve1$s
## Send it back to the count space


