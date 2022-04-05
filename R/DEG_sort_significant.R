library(dplyr)

rm(list=ls())

# method = 'vanilla'
method = 'change'
# sort_method = 'bayes'
sort_method = 'fdr'
marker_genes <- F
# marker_genes <- T
no_CellTypes <- F

dir.home <- getwd()
dir.data <- paste(dir.home, '/results/DEG_analysis', sep = '')
samp <- list.files(dir.data)
samp <- samp[1] # change this if calculate for a different sample
if (marker_genes == T){
  samp <- paste(samp, '/marker_genes', sep = '')
} else if (no_CellTypes == T){
  samp <- paste(samp, '/no_CellTyping', sep = '')
}
if (method == 'vanilla'){
  dir.data <- paste(dir.data, samp, 'vanilla_mode', sep = '/')
} else if (method == 'change'){
  dir.data <- paste(dir.data, samp, 'change_mode', sep = '/')
}
# dir.data <- paste(dir.data, '/change_delta_001/FCs', sep = '')
dir.out <- paste(dir.data, '/', sort_method, '_sorted', sep = '')

if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print('dir.out was created')
}

lists <- list.files(path = dir.data, full.names = T, pattern = 'DEGs_')

X <- c(1:length(lists))
for (i in 1:length(lists)){
  X[i] <- lapply(lists[i], read.table, header = T, row.names = 1, sep = ',')
}

# Sort out the significant genes
for (i in 1:length(lists)){
  if (method == 'change'){
    # X[[i]] <- X[[i]][which(X[[i]]$proba_not_de <= 0.05),]
    if (sort_method == 'fdr'){
      X[[i]] <- X[[i]][which(X[[i]]$is_de_fdr_0.05 == 'True'),]
    } else if (sort_method == 'bayes'){
      X[[i]] <- X[[i]][which(abs(X[[i]]$bayes_factor) >= 3),]
    }
  } else if (method == 'vanilla'){
    # X[[i]] <- X[[i]][which(X[[i]]$pval <= 0.05),]
    X[[i]] <- X[[i]][which(abs(X[[i]]$bayes_factor) >= 3),]
    # X[[i]] <- X[[i]][which(abs(X[[i]]$bayes_factor) >= 2.3),]
    # X[[i]] <- X[[i]][which(X[[i]]$proba_m1 <= 0.05 | X[[i]]$proba_m2 <= 0.05),]
  }
}


# Calculate FC
if (method == "vanilla"){
  for (i in 1:length(lists)){
    X[[i]]$'logFC' <- log2(X[[i]]$scale1/X[[i]]$scale2)
  }
}

# subset columns of interest
cnms <- colnames(X[[1]]) # colnames of all dataframes should be the same
for (i in 1:length(lists)){
  if (method == 'change'){
    X[[i]] <- X[[i]][,which(cnms %in% c('proba_de', 'proba_not_de', 'bayes_factor', 'lfc_mean', 'is_de_fdr_0.05'))]
  } else if (method == 'vanilla'){
    X[[i]] <- X[[i]][,which(cnms %in% c('proba_m1', 'proba_m2', 'bayes_factor', 'logFC'))]
  }
}


# Write to out
outnames <- paste('significant_', sapply(strsplit(lists, '/'), tail, 1), sep = '')
for (i in 1:length(lists)){
  if (length(rownames(X[[i]])) < 1){
    next
  }
  write.csv(X[i], paste(dir.out, outnames[i], sep = '/'))
}

#########################
### Create combined table of DEGs per cell type for NicheNet analysis
#########################
list_DEGs_NicheNet <- function(X, lists){
  groups <- sapply(strsplit(lists, '/'), tail, 1)
  groups <- sapply(strsplit(groups, '\\.'), '[[', 1)
  groups <- paste(sapply(strsplit(groups, '_'), '[[', 5),
                  sapply(strsplit(groups, '_'), '[[', 6), sep = '_')
  
  if (length(groups) != length(unique(sort(groups)))){
    stop('The groups of interactions are not unique')
  }
  
  DEGs_maxlength <- c()
  groups_x <- c()
  for (i in 1:length(X)){
    if (length(rownames(X[[i]])) < 1){
      next
    }
    DEGs_maxlength <- c(DEGs_maxlength, length(rownames(X[[i]])))
    groups_x <- c(groups_x, groups[i])
  }
  DEGs_Ngroups <- length(DEGs_maxlength)
  DEGs_maxlength <- max(DEGs_maxlength)
  DEGs <- matrix(NA, nrow = DEGs_maxlength, ncol = DEGs_Ngroups)
  # remove(DEGs_maxlength, DEGs_Ngroups)
  z <- 1
  for (i in 1:length(X)){
    if (length(rownames(X[[i]])) < 1){
      next
    }
    # colnames(DEGs)[z] <- groups[i]
    DEGs[1:length(rownames(X[[i]])),z] <- rownames(X[[i]])
    z <- z+1
  }
  colnames(DEGs) <- groups_x
  return(DEGs)
}

DEGs <- list_DEGs_NicheNet(X, lists)
write.csv(DEGs, paste(dir.out, '/DEGs_Sick_vs_Healthy.csv', sep = ''), row.names = F)



