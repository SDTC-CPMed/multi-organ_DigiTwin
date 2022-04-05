library(Seurat)
library(dplyr)
library(Matrix)
library(R.utils)
library(patchwork)


rm(list=ls())
dir.home <- getwd()
dir.home <- paste(dir.home, '/DGE_data/', sep = '')
dir.data <- dir.home
dir.out <- paste(dir.home, 'sorted_DGEs/', sep = '')
if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print('dir.out was created')
}

projectname <- sapply(strsplit(dir.home, '/'), '[[', 4)
files <- list.files(dir.data, pattern = 'txt.gz')
tissues <- sapply(strsplit(files, '_'), '[[', 1)

# tissue <- tissues[1]
for (tissue in tissues){
  X <- read.table(paste(dir.data, list.files(dir.data, pattern = tissue), sep = ''), header = T, sep = '\t', row.names = 1)
  X[1:5,1:5]
  X_backup <- X
  
  ### Seurat QC and data sorting
  X.seu <- CreateSeuratObject(X, min.cells = ceiling(length(colnames(X))*0.01))
  X.seu[["percent.mito"]] <- PercentageFeatureSet(X.seu, pattern = "^mt-")
  X.seu <- subset(X.seu, subset = nFeature_RNA > 200 & nCount_RNA > 400 & percent.mito < 20)

  # VlnPlot(X.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  # plot1 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "percent.mito")
  # plot2 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # plot1 + plot2
  
  ### Sort the data based on results from previous graphs
  if (tissue == 'Joint'){
    print(tissue)
  } else if (tissue == 'Lung'){
    print(tissue)
    X.seu <- subset(X.seu, subset = nCount_RNA < 2000)
  } else if (tissue %in% c('Spleen')){
    print(tissue)
    X.seu <- subset(X.seu, subset = nCount_RNA < 6000)
  } else if (tissue %in% c('Muscle', 'Skin')){
    print(tissue)
    X.seu <- subset(X.seu, subset = nCount_RNA < 7000)
  } else {
    print('tissue specific sorting conditions not defined')
  }

  ### Extract expression matrix
  X <- as.data.frame(as.matrix(GetAssayData(X.seu, slot = "counts")))
  X[1:5,1:5]
  write.csv(X, paste(dir.out,  paste(tissue, '_sorted_expression_matrix.csv', sep = ''), sep = ''), quote = F)

}

