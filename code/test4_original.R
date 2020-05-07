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
low_to_high2 <- function(dat, pca, replaced_x) {
  mu <- colMeans(dat)
  
  ncomponents <- 15
  dat_hat <- replaced_x[,1:ncomponents] %*% t(pca$rotation[,1:ncomponents])
  dat_hat <- scale(dat_hat, center = -mu, scale = FALSE)
  
  return(dat_hat)
}

pca <- prcomp(t(logcounts(se)), center = TRUE, scale = FALSE, rank. = 15)

new_count_matrix <- low_to_high2(t(logcounts(se)), pca, sds@curves$curve1$s)

# Wrong way to check
# check <- calculatePCA(t(new_count_matrix), ncomponents = 15)
# hist(check - rd, main = "sds curve to count matrix and back ")


new_count_matrix <- low_to_high2(t(logcounts(se)), pca, pca$x[, 1:15])
check <- prcomp(new_count_matrix, center = TRUE, scale = FALSE, rank. = 15)
hist(abs(check$x[, 1:15]) - abs(pca$x[, 1:15]), main = "Count matrix to pca to count matrix to pca")

for(i in 1:15) plot(pca$x[,i], check$x[,i])

# FitGAM
sce <- fitGAM(counts = se)
saveRDS(sce, file = "fitGAM_21.rds")
GAM <- read_rds("fitGAM_21.rds")

startRes <- startVsEndTest(GAM)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[1]]
plotSmoothers(sce, counts, gene = sigGeneStart)

predictions <- predictSmooth(GAM, gene = oStart[1:20], nPoints = 40) 