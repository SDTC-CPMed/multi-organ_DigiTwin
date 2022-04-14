# Sandra Lilja
# R version 4.0.4
#'
#' Subset the lists of DEGs only to include significant DEGs. 
#' If method == 'change', the significant DEGs are those with is_de_fdr_0.05 == 'True'
#' If method == 'vanilla', the significant DEGs are those with abs(bayes_factor) >= 3
#' If method == 'vanilla', the fold changes will also be calculated by log2(scale1/scale2). 
#' 
#' @param DEG_files The filenames (including paths) of the csv input file containing all DEG information
#' @param outdir_DEG The path to the input files, used for creation of output path.
#' @param method Optional parameter defining the scVI method based on which DEGs were calculated. Can be either 'change' or 'vanilla'. Default: methods = 'change'. 
#' 
#' @export
#'

library(dplyr)

DEG_sort_significant <- function(DEG_files, outdir_DEG, method = 'change'){
  dir.out <- paste(outdir_DEG, '/', 'subset_significant', sep = '')
  if (dir.exists(dir.out)==FALSE){
    dir.create(dir.out)
  }
  
  X <- c(1:length(DEG_files))
  for (i in 1:length(DEG_files)){
    X[i] <- lapply(DEG_files[i], read.table, header = T, row.names = 1, sep = ',')
  }
  
  # Sort out the significant genes
  for (i in 1:length(DEG_files)){
    if (method == 'change'){
      X[[i]] <- X[[i]][which(X[[i]]$is_de_fdr_0.05 == 'True'),]
    } else if (method == 'vanilla'){
      X[[i]] <- X[[i]][which(abs(X[[i]]$bayes_factor) >= 3),]
      # Calculate FC (if method = 'vanilla')
      X[[i]]$'logFC' <- log2(X[[i]]$scale1/X[[i]]$scale2)
    }
  }
  
  # subset columns of interest
  cnms <- colnames(X[[1]]) # colnames of all dataframes should be the same
  for (i in 1:length(DEG_files)){
    if (method == 'change'){
      X[[i]] <- X[[i]][,which(cnms %in% c('proba_de', 'proba_not_de', 'bayes_factor', 'lfc_mean', 'is_de_fdr_0.05'))]
    } else if (method == 'vanilla'){
      X[[i]] <- X[[i]][,which(cnms %in% c('proba_m1', 'proba_m2', 'bayes_factor', 'logFC'))]
    }
  }
  
  # Write to out
  outnames <- paste('significant_', sapply(strsplit(DEG_files, '/'), tail, 1), sep = '')
  for (i in 1:length(DEG_files)){
    if (length(rownames(X[[i]])) < 1){
      next
    }
    write.csv(X[i], paste(dir.out, outnames[i], sep = '/'))
  }
  
  # #########################
  # ### Create combined table of DEGs per cell type for NicheNet analysis
  # #########################
  # list_DEGs_NicheNet <- function(X, DEG_files){
  #   groups <- sapply(strsplit(DEG_files, '/'), tail, 1)
  #   groups <- sapply(strsplit(groups, '\\.'), '[[', 1)
  #   groups <- paste(sapply(strsplit(groups, '_'), '[[', 5),
  #                   sapply(strsplit(groups, '_'), '[[', 6), sep = '_')
  #   
  #   if (length(groups) != length(unique(sort(groups)))){
  #     stop('The groups of interactions are not unique')
  #   }
  #   
  #   DEGs_maxlength <- c()
  #   groups_x <- c()
  #   for (i in 1:length(X)){
  #     if (length(rownames(X[[i]])) < 1){
  #       next
  #     }
  #     DEGs_maxlength <- c(DEGs_maxlength, length(rownames(X[[i]])))
  #     groups_x <- c(groups_x, groups[i])
  #   }
  #   DEGs_Ngroups <- length(DEGs_maxlength)
  #   DEGs_maxlength <- max(DEGs_maxlength)
  #   DEGs <- matrix(NA, nrow = DEGs_maxlength, ncol = DEGs_Ngroups)
  #   # remove(DEGs_maxlength, DEGs_Ngroups)
  #   z <- 1
  #   for (i in 1:length(X)){
  #     if (length(rownames(X[[i]])) < 1){
  #       next
  #     }
  #     # colnames(DEGs)[z] <- groups[i]
  #     DEGs[1:length(rownames(X[[i]])),z] <- rownames(X[[i]])
  #     z <- z+1
  #   }
  #   colnames(DEGs) <- groups_x
  #   return(DEGs)
  # }
  # 
  # DEGs <- list_DEGs_NicheNet(X, DEG_files)
  # write.csv(DEGs, paste(dir.out, '/DEGs_Sick_vs_Healthy.csv', sep = ''), row.names = F)
  
}




