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

species <- "mouse"
projectname <- sapply(strsplit(dir.home, '/'), '[[', 4)
files <- list.files(dir.data, pattern = 'txt.gz')
tissues <- sapply(strsplit(files, '_'), '[[', 1)
# tissues <- 'allTissues'

# tissue <- tissues[1]
for (tissue in tissues){
  X <- read.table(paste(dir.data, list.files(dir.data, pattern = tissue), sep = ''), header = T, sep = '\t', row.names = 1)
  X[1:5,1:5]
  X_backup <- X
  
  ### Seurat QC and data sorting
  ### if needed, here I can also change the minimum requirement of number of reads
  X.seu <- CreateSeuratObject(X, min.cells = ceiling(length(colnames(X))*0.01))
  if (species == "mouse"){
    X.seu[["percent.mito"]] <- PercentageFeatureSet(X.seu, pattern = "^mt-")
    X.seu <- subset(X.seu, subset = nFeature_RNA > 200 & nCount_RNA > 400 & percent.mito < 20)
  } else if(species == "human"){
    X.seu[["percent.mito"]] <- PercentageFeatureSet(X.seu, pattern = "^MT-")
    X.seu <- subset(X.seu, subset = nFeature_RNA > 200 & nCount_RNA > 400 & percent.mito < 20)
  } else {
    print("Which species?")
  }
  VlnPlot(X.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  
  plot1 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "percent.mito")
  plot2 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  vv <- X.seu@meta.data$nFeature_RNA/X.seu@meta.data$nCount_RNA
  plot(density(vv))
  
  head(X.seu@meta.data)
  min(X.seu@meta.data[,3])
  max(X.seu@meta.data[,3])
  min(X.seu@meta.data[,2])
  max(X.seu@meta.data[,4])
  length(which(X.seu@meta.data[,2]>=6000)==TRUE)
  length(which(X.seu@meta.data[,3]>=2500)==TRUE)
  length(rownames(X.seu@meta.data))
  
  ## Removing outliers based on quantile would remove too much!
  # qnt <- quantile(X.seu@meta.data[,2], probs = c(0.25, 0.75))
  # H <- 1.5 * IQR(X.seu@meta.data[,2])
  # qnt[1]-H
  # qnt[2]+H
  
  ### Sort the data based on results from previous graph
  ### Note, the readme file sais that no outliers were removed fro tissue == 'HC'
  if (tissue == 'Joint'){
    print(tissue)
  } else if (tissue == 'Lung'){
    print(tissue)
    X.seu <- subset(X.seu, subset = nCount_RNA < 2000)
  } else if (tissue %in% c('Spleen', 'allTissues')){
    print(tissue)
    X.seu <- subset(X.seu, subset = nCount_RNA < 6000)
  } else if (tissue %in% c('Muscle', 'Skin')){
    print(tissue)
    X.seu <- subset(X.seu, subset = nCount_RNA < 7000)
  } else {
    print('tissue specific sorting conditions not defined')
  }

  ### Extract the Seurat results
  X <- as.data.frame(as.matrix(GetAssayData(X.seu, slot = "counts")))
  X[1:5,1:5]
  # write.table(X, paste(dir.out, paste(tissue, '_sorted_expression_matrix.txt', sep = ''), sep = '/'),
  #             row.names = T, col.names = T, sep = '\t', quote = F)
  write.csv(X, paste(dir.out,  paste(tissue, '_sorted_expression_matrix.csv', sep = ''), sep = ''), quote = F)
  
  
  ### QC post-sorting
  if (length(grep('qc_results.csv', list.files(dir.home))) == 1){
    QC <- as.matrix(read.csv(paste(dir.home, 'qc_results.csv', sep = '/'), sep = ',', header = T))
    QC <- rbind(QC, NA)
  } else {
    QC <- matrix(data = NA, ncol = 11, nrow = 1)
    colnames(QC) <- c('project.name', 'date', 'tissue', 'Ncells', 'Ncells above 500 genes', 'Ncells above 1000 genes', 'mean Ngenes per cell', 'min Ngenes per cell', 'max Ngenes per cell', 'max NUMI per cell', 'mean NUMI per cell')
  }

  # head(QC)
  # head(X.seu@meta.data)
  QC[nrow(QC),1] <- projectname
  QC[nrow(QC),2] <- as.character(Sys.Date())
  QC[nrow(QC),3] <- tissue
  QC[nrow(QC),4] <- length(X.seu@meta.data[,1])
  QC[nrow(QC),5] <- length(which(X.seu@meta.data[,3]>500)==TRUE)
  QC[nrow(QC),6] <- length(which(X.seu@meta.data[,3]>1000)==TRUE)
  QC[nrow(QC),7] <- mean(X.seu@meta.data[,3]) # meanGene
  QC[nrow(QC),8] <- min(X.seu@meta.data[,3]) # minGene
  QC[nrow(QC),9] <- max(X.seu@meta.data[,3]) # maxGene
  QC[nrow(QC),10] <- max(X.seu@meta.data[,2]) # maxUMI
  QC[nrow(QC),11] <- mean(X.seu@meta.data[,2]) # meanUMI

  write.csv(QC, paste(dir.home, 'qc_results.csv', sep = '/'), row.names = F)
  
}

# filenames <- list.files(dir.out, pattern = '.txt', full.names = T)
# for (filename in filenames){
#   gzip(filename)
# }

