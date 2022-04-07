# Sandra Lilja
# R version 4.0.4
#'
#' Check the quality of the data and remove outliers
#' 
#' @param filename The filename of the csv.gz input file with cells in columns and genes in rows
#' @param organ The organ as a string (eg., "Joint"). Based on this parameter, outliers will be romoved as manually defined for our data. 
#' 
#' @export
#'

library(Seurat)

sc_data_quality_sorting <- function(filename, organ){
  # X <- read.table(paste(dir.data, list.files(dir.data, pattern = organ), sep = ''), header = T, sep = '\t', row.names = 1)
  X <- read.table(filename, header = T, sep = '\t', row.names = 1)
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
  if (organ == 'Joint'){
  } else if (organ == 'Lung'){
    X.seu <- subset(X.seu, subset = nCount_RNA < 2000)
  } else if (organ %in% c('Spleen')){
    X.seu <- subset(X.seu, subset = nCount_RNA < 6000)
  } else if (organ %in% c('Muscle', 'Skin')){
    X.seu <- subset(X.seu, subset = nCount_RNA < 7000)
  } else {
    print('organ specific sorting conditions not defined')
  }
  
  ### Extract expression matrix
  X <- as.data.frame(as.matrix(GetAssayData(X.seu, slot = "counts")))
  X[1:5,1:5]
  # write.csv(X, paste(dir.out,  paste(organ, '_sorted_expression_matrix.csv', sep = ''), sep = ''), quote = F)
  write.csv(X, paste('data/sorted_DGEs/', organ, '_sorted_expression_matrix.csv', sep = ''), quote = F)
  
}

