library(SingleCellExperiment)
library(scater)
library(slingshot)
library(tradeSeq)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(remotes)
library(readr)
library(stats)
load("../data/se_paul.rda")

set.seed(21)
counts(se) <- as.matrix(counts(se))
se <- logNormCounts(se)
se <- runUMAP(se)
se <- runPCA(se, ncomponents = 15)
rd <- reducedDim(se, "PCA")
cl <- kmeans(rd, 7)$cluster
se$kmeans <- cl

se <- slingshot(se, clusterLabels = "kmeans", reducedDim = "PCA", start.clus = 3)
sds <- SlingshotDataSet(se)

# you don't need to define pca as input to the function and inside the functio 

low_to_high <- function(dat, pca, replaced_x, ncomponents) {
  
  dat_hat <- replaced_x[,1:ncomponents] %*% t(pca$rotation[,1:ncomponents])
  dat_hat <- scale(dat_hat, center = TRUE, scale = FALSE)
  
  return(dat_hat)
}

pca <- prcomp(t(logcounts(se)), center = TRUE, scale = FALSE, rank. = 15)
new_count_matrix <- low_to_high(t(logcounts(se)), pca, sds@curves$curve1$s, 15)

# Wrong way to check
# check <- calculatePCA(t(new_count_matrix), ncomponents = 15)
# hist(check - rd, main = "sds curve to count matrix and back ")


new_count_check <- low_to_high(t(logcounts(se)), pca, pca$x[, 1:15])
check <- prcomp(new_count_check, center = TRUE, scale = FALSE, rank. = 15)
hist(abs(check$x[, 1:15]) - abs(pca$x[, 1:15]), main = "Count matrix to pca to count matrix to pca")

for(i in 1:15) plot(pca$x[,i], check$x[,i])

# FitGAM
# sce <- fitGAM(counts = se)
# saveRDS(sce, file = "fitGAM_21.rds")
GAM <- read_rds("fitGAM_21.rds")

startRes <- startVsEndTest(GAM)
oStart <- order(startRes$waldStat, decreasing = TRUE)

for(i in 1:9){
  sigGeneStart <- names(GAM)[oStart[1]]
  plotSmoothers(GAM, counts(GAM), gene = sigGeneStart)
}

par(mfrow = c(3,3))
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[1]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[2]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[3]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[4]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[5]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[6]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[7]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[8]])
plotSmoothers(GAM, t(new_count_matrix), gene = names(GAM)[oStart[9]])

predictions <- predictSmooth(GAM, gene = oStart[1:20], nPoints = 40) 