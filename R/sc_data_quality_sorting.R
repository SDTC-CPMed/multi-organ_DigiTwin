# Sandra Lilja
# R version 4.0.4
#'
#' Check the quality of the data and remove outliers
#' 
#' @param filename The filename of the csv.gz input file with cells in columns and genes in rows
#' 
#' @export
#'

library(Seurat)

sc_data_quality_sorting <- function(filenames){
  X <- read.table(filename, header = T, sep = '\t', row.names = 1)
  
  ### Seurat QC and data sorting
  X.seu <- CreateSeuratObject(X, min.cells = ceiling(length(colnames(X))*0.01))
  X.seu[["percent.mito"]] <- PercentageFeatureSet(X.seu, pattern = "^mt-")
  X.seu <- subset(X.seu, subset = nFeature_RNA > 200 & nCount_RNA > 400 & percent.mito < 20)
  
  # VlnPlot(X.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  # plot1 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "percent.mito")
  # plot2 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # plot1 + plot2
  
  ### Sort the data based on results from previous graphs
  X.seu <- subset(X.seu, subset = nCount_RNA < 6000)

  ### Extract expression matrix
  X <- as.data.frame(as.matrix(GetAssayData(X.seu, slot = "counts")))
  write.csv(X, paste('data/sorted_DGEs/sorted_expression_matrix.csv', sep = ''), quote = F)
  
}

