library(Seurat)
library(dplyr)
library(Matrix)
library(plyr)
library(R.utils)

rm(list=ls())
dir.home <- getwd()
dir.data <- paste(dir.home, '/scData/', sep = '')
dir.out <- paste(dir.home, '/DGE_data/', sep = '')
if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print('dir.out was created')
}

species <- "mouse"

dirs <- list.files(path = dir.data, full.names = T, pattern = 'NS500340|-')
dirs <- list.files(path = dirs, full.names = T, pattern = 'CIA')


X <- c(1:length(dirs))
for (i in 1:length(dirs)){
  X[i] <- lapply(paste(dirs[i], '/umi.dge.txt.gz', sep = ''), read.table, header = T, row.names = 1, sep = '\t')
}
X[[1]][1:5,1:5]
X[[6]][1:5,1:5]
X[[39]][1:5,1:5]

sampnames <- sapply(strsplit(dirs, '/'), tail, 1)
## restructure names of old blood samples
sampnames[9]
sampnames[1:8]
sampnames[1:8] <- gsub('_mouse_', '_', sampnames[1:8])
sampnames[1:8] <- gsub('_healthy_', '_H', sampnames[1:8])
sampnames[1:8] <- gsub('_sick_', '_S', sampnames[1:8])
sampnames[1:8] <- paste(sapply(strsplit(sampnames[1:8], '_'), '[[', 3),
      sapply(strsplit(sampnames[1:8], '_'), '[[', 2),
      sapply(strsplit(sampnames[1:8], '_'), '[[', 1),
      sapply(strsplit(sampnames[1:8], '_'), '[[', 4),
      sapply(strsplit(sampnames[1:8], '_'), '[[', 5), sep = '_')
##
sampnames <- paste(sapply(strsplit(sampnames, '_'), '[[', 3),
                   sapply(strsplit(sampnames, '_'), '[[', 1),
                   sapply(strsplit(sampnames, '_'), '[[', 2), sep = '_')
length(sampnames) == length(unique(sort(sampnames)))
length(sampnames) == length(X)


for (i in 1:length(sampnames)){
  colnames(X[[i]]) <- paste(sampnames[i], colnames(X[[i]]), sep = '_')
}
X[[1]][1:5,1:5]
X[[5]][1:5,1:5]
X[[9]][1:5,1:5]


as.data.frame(lapply(X, dim))

X_backup <- X
# X <- X_backup
# for (i in 1:length(X)){
#   print(i)
#   print(grep('_', rownames(X[[i]])))
# }
# rownames(X[[42]])[13715]
# grep('RMST', rownames(X[[42]]))

if (species == "mouse"){
  for (i in 1:length(X)){
    X.seu <- CreateSeuratObject(X[[i]], min.cells = 3)
    X.seu[["percent.mito"]] <- PercentageFeatureSet(X.seu, pattern = "^mt-")
    X.seu <- subset(X.seu, subset = nFeature_RNA > 200 & nCount_RNA > 400 & percent.mito < 20)
    X[[i]] <- X[[i]][,c(rownames(X.seu@meta.data))]
  }
} else if(species == "human"){
  for (i in 1:length(X)){
    X.seu <- CreateSeuratObject(X[[i]], min.cells = 3)
    X.seu[["percent.mito"]] <- PercentageFeatureSet(X.seu, pattern = "^MT-")
    X.seu <- subset(X.seu, subset = nFeature_RNA > 200 & nCount_RNA > 400 & percent.mito < 20)
    X[[i]] <- X[[i]][,c(rownames(X.seu@meta.data))]
  }
} else {
  print("Which species?")
}
X[[5]][1:5,1:5]
as.data.frame(lapply(X, dim))  

head(colnames(X[[1]]))
head(rownames(X[[1]]))
ncells <- sum(as.data.frame(lapply(X, dim))[2,])  

for (i in 1:length(X)){
  X[[i]]$rn <- rownames(X[[i]])
}
X[[5]][1:5,(length(X[[5]])-5):length(X[[5]])]


## Separate the different tissues into different DGE-files
# tissues <- unique(sort(sapply(strsplit(sampnames, '_'), '[[', 1)))
tissues <- 'allTissues'
# tissue = tissues[1]
for (tissue in tissues){
  if (tissue == 'allTissues'){
    dfs <- 1:length(sampnames)
    Y <- X
  } else if (tissue != 'allTissues'){
    dfs <- grep(tissue, sampnames)
    Y <- X[c(dfs)]
  }
  ncells_tissue <- sum(as.data.frame(lapply(Y, dim))[2,])  
  Y <- join_all(Y, by = 'rn', type = 'full')
  rownames(Y) <- Y$rn
  Y[1:5,1:5]
  Y$rn <- NULL
  if (tissue != 'allTissues'){
    if ((unique(sort(sapply(strsplit(colnames(Y), '_'), '[[', 1))) == tissue) == FALSE){
      print('Error: Tissue sorting need revising')
      break
    }
  }
  if ((length(colnames(Y))+length(dfs) == ncells_tissue)==FALSE){
    print('Error: Not all cells were joined')
    break
  }
  if (any(is.na(Y)) == TRUE){
    Y[is.na(Y)] <- 0
  } else {
    print(any(is.na(Y)))
  }
  write.table(Y, paste(dir.out, tissue, '_UMI_expression_matrix.txt', sep = ''), col.names = T, row.names = T, sep = '\t', quote = F)
}


filenames <- list.files(dir.out, pattern = 'UMI_expression_matrix.txt', full.names = T)
filenames <- subset(filenames, !grepl('txt.gz', filenames))
for (filename in filenames){
  gzip(filename)
}
