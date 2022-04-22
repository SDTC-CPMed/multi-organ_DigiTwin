# Sandra Lilja
# R version 4.0.4
#'
#' Plot the expression level and/or FC of a set of genes 
#' 
#' @param exprdata The normalized single-cell expression matrix
#' @param clusts The clusters and cell type information, containing the column 'CellType'
#' @param orth A list of mouse genes and their human orthologs
#' @param DEGs  a matrix containing lists of DEGs for each cell type and organ
#' @param dir.out  output directory
#' @param ligand_target
#' @param lr_network
#' 
#' @export
#'

library(nichenetr)
library(dplyr)

# ct <- 'cluster_ids.csv' # define cell type translation file
# orth <- 'orthologous_translation_file.txt' # define orthologus translation file
# 
# dir.data <- '../data/allTissues_tissue-sample-BatchRemoval'
# dir.ct <- '../data/clusters_final_out'
# dir.degs <- '../data/DEG_analysis/fdr_sorted'
# dir.orth <- 'data'
# dir.out <- '../data/NicheNet_analysis'
# if (dir.exists(dir.out)==FALSE){
#   dir.create(dir.out, recursive = T)
#   print('dir.out was created')
# }
# 
# # scVI adjusted scRNA-seq mouse data 
# exprdata <- read.csv(paste(dir.data, "/normalized_expression_matrix.csv", sep = ''), row.names = 1)
# # orthologous translation file
# orth <- read.table(paste(dir.orth, orth, sep = '/'), header = T)


list_orthologs  <- function(mouse_gene_list, orth){
  orth_sub = orth[which(orth$Mouse.gene.name %in% mouse_gene_list == TRUE),]
  # human_gene_orthologs <- as.character(orth$Gene.name[which(orth$Mouse.gene.name %in% mouse_gene_list == TRUE)])
  return(orth_sub)
}
find_orthologs_mouse_to_human  <- function(mouse_gene_list, orth){
  human_gene_orthologs <- as.character(orth$Gene.name[which(orth$Mouse.gene.name %in% mouse_gene_list == TRUE)])
  return(human_gene_orthologs)
}

# NicheNet analysis steps referred to below are from this vignette:
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md

# STEP 1: define expressed genes in sender and receiver cell population
# s = cell IDs sender
# t = cell IDs target
define_expressed_genes <- function(exprdata, s = colnames(exprdata), t = colnames(exprdata), exp_cutoff = 0.1){
  
  # exp_cutoff <- 0.0001
  
  # reverse log10 tranformation
  ### Why ?
  expressed_genes_sender <- 10**exprdata[,colnames(exprdata) %in% s]-1
  expressed_genes_target <- 10**exprdata[,colnames(exprdata) %in% t]-1
  # expressed_genes_sender <- exprdata[,colnames(exprdata) %in% s]
  # expressed_genes_target <- exprdata[,colnames(exprdata) %in% t]
  # apply NicheNet recommended expression threshold
  temp <- log2(rowMeans(expressed_genes_sender)+1) >= exp_cutoff
  expressed_genes_sender <- rownames(expressed_genes_sender)[temp]
  temp <- log2(rowMeans(expressed_genes_target)+1) >= exp_cutoff
  expressed_genes_target <- rownames(expressed_genes_target)[temp]
  print(paste("n genes expressed in sender:", length(expressed_genes_sender), " n genes expressed in target:", length(expressed_genes_target), sep=" "))
  
  exp_genes <- cbind(c(expressed_genes_sender, rep(NA, times = nrow(exprdata)-length(expressed_genes_sender))),
                     c(expressed_genes_target, rep(NA, times = nrow(exprdata)-length(expressed_genes_target))))
  colnames(exp_genes) <- c("Sender", "Target")
  
  return(exp_genes)
}

# STEP 2: define gene set of interest and a background of genes (receiver)
# STEP 3: define a set of potential ligands (sender)
# exp_genes = derived from step 1
genes_of_interest <- function(sender_degs, target_degs, exp_genes){
  sender_interest_gene <- sender_degs[sender_degs %in% colnames(ligand_target)]
  sender_interest_gene <- sender_interest_gene[sender_interest_gene %in% lr_network[lr_network[,2] %in% exp_genes[,2],1]] # only select ligands that target a receptor expressed in receiver
  target_interest_gene <- target_degs[target_degs %in% rownames(ligand_target)]
  
  print(paste("n sender DEGs not in expressed genes:", sum(!(sender_interest_gene %in% exp_genes[,1])), sep=" "))
  sender_interest_gene <- sender_interest_gene[sender_interest_gene%in%exp_genes[,1]]
  print(paste("n target DEGs not in expressed genes:", sum(!(target_interest_gene %in% exp_genes[,2])), sep=" "))
  target_interest_gene <- target_interest_gene[target_interest_gene%in%exp_genes[,2]]
  
  int_genes <- cbind(c(sender_interest_gene, rep(NA, times = nrow(DEGs)-length(sender_interest_gene))),
                     c(target_interest_gene, rep(NA, times = nrow(DEGs)-length(target_interest_gene))))
  colnames(int_genes) <- c("Sender", "Target")
  return(int_genes)
}

# STEP 4: perform NicheNets ligand activity analysis on the gene set of interest
nichenet_analysis <- function(exp_genes, int_genes, ligand_target){
  
  ligand_activity2 <- predict_ligand_activities2(geneset = int_genes[!is.na(int_genes[,2]),2],
                                                 background_expressed_genes = exp_genes[!is.na(exp_genes[,2]),2],
                                                 ligand_target_matrix = ligand_target,
                                                 potential_ligands = unique(int_genes[!is.na(int_genes[,1]),1]),
                                                 single = T)
  ligand_activity <- predict_ligand_activities(geneset = int_genes[!is.na(int_genes[,2]),2],
                                               background_expressed_genes = exp_genes[!is.na(exp_genes[,2]),2],
                                               ligand_target_matrix = ligand_target,
                                               potential_ligands = unique(int_genes[!is.na(int_genes[,1]),1]),
                                               single = T)
  # rank ligands based on ligand activity (pearson correlation coefficient)
  # comment below if predict_ligand_activities2
  ligand_activity <- ligand_activity %>% arrange(-pearson)
  # ligand_activity <- as.matrix(ligand_activity)
  # ligand_activity <- ligand_activity[order(ligand_activity[,4], decreasing = T),]
  # end comment
  return(list(ligand_activity, ligand_activity2))
}

extractCorrInfo <- function(setting, ligand_target_matrix){
  #browser()
  PredictionG=setting$response
  TargetG=ligand_target_matrix[,setting$from]
  names(TargetG)<-rownames(ligand_target_matrix)
  
  response_df = tibble(gene = names(TargetG), response =TargetG)
  prediction_df = tibble(gene = names(PredictionG), prediction = PredictionG)
  combined = inner_join(response_df, prediction_df, by = "gene")
  corp=cor(combined$prediction, combined$response)
  
  return(list(from=setting$from, info=combined, cor=corp ))
  
}

# "ranking" parameter contains rankings per couple of cell types
# info parameter contains infomation needed to aggregate information over several cell types

predict_ligand_activities2 <- function (geneset, background_expressed_genes,ligand_target_matrix, potential_ligands, single = TRUE,...){
  
  setting = list(geneset) %>% lapply(convert_gene_list_settings_evaluation, 
                                     name = "gene set", ligands_oi = potential_ligands, background = background_expressed_genes)
  if (single == TRUE) {
    settings_ligand_prediction = setting %>% convert_settings_ligand_prediction(all_ligands = potential_ligands, 
                                                                                validation = FALSE, single = TRUE)
    
    # browser()
    ligand_importance_Info = settings_ligand_prediction %>% lapply(extractCorrInfo, 
                                                                   ligand_target_matrix = ligand_target_matrix)
    
    names(ligand_importance_Info)<-potential_ligands
    
    results_ranking=ligand_importance_Info%>%lapply(function(x) data.frame(from=as.character(x$from), cor=x$cor) )%>%
      bind_rows()%>%arrange(desc(cor))
    
    
  }
  
  return(list(info=ligand_importance_Info, ranking=results_ranking))
  
}


NicheNet_analysis_main <- function(exprdata, clusts, orth, DEGs, dir.out, ligand_target, lr_network){
  # load the expr data
  exprdata <- t(exprdata) # as NicheNet requires the transposed input
  # exprdata[1:5,1:3]
  # Cell types
  cell_types <- as.matrix(clusts)
  # remove any unused column. CellType_tissue must be in the fifth column
  CellType_tissue <- as.vector(sub(' ', '_', paste(gsub(' ', '-', cell_types[,'CellType_markers']), cell_types[,'tissue'])))
  cell_types <- cbind(cell_types, CellType_tissue)
  remove(CellType_tissue)
  # Get DEGs for cell types
  DEGs <- as.matrix(DEGs)
  colnames(DEGs) <- gsub('\\.', '-', colnames(DEGs))
  
  # # Get ligand-target matrix
  # ligand_target <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds")) # targets = rows, ligands = columns
  # 
  # # Get ligand-receptor interactions
  # lr_network = as.matrix(readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds")))
  
  # CONVERT GENES TO HUMAN ORTHOLOGS
  #########################################################################################################
  # Ensure that the orthologous file contains distinct human/mouse translations
  orth <- orth[,c('Gene.name', 'Mouse.gene.name')]
  orth <- distinct(orth)
  # for each list of DEGs, translate to human orthologs
  for (col in 1:length(colnames(DEGs))){
    # col <- 8
    mouse_gene_list <- DEGs[,col]
    orth_x <- list_orthologs(mouse_gene_list, orth)
    if (exists('orth_sub') == FALSE){
      if (length(rownames(orth_x)) != 0){
        orth_sub <- orth_x
      }
    } else {
      if (length(rownames(orth_x)) != 0){
        orth_sub <- rbind(orth_sub, list_orthologs(mouse_gene_list, orth))
      }
    }
    human_gene_orthologs <- unique(sort(find_orthologs_mouse_to_human(mouse_gene_list, orth)))
    DEGs[,col] <- NA
    if (length(human_gene_orthologs) != 0){
      DEGs[1:length(human_gene_orthologs),col] <- human_gene_orthologs
    }
  }
  # remove samples where no orthologs genes could be found
  coln <- c()
  for (i in 1:length(colnames(DEGs))){
    if (length(unique(DEGs[,i])) == 1){
      if (is.na(unique(DEGs[,i]))){
        next
      }
    }
    else {
      coln <- c(coln, colnames(DEGs)[i])
      if (exists('degs_x')){
        degs_x <- cbind(degs_x, DEGs[,i])
      }
      else {
        degs_x <- DEGs[,i]
      }
    }
  }
  colnames(degs_x) <- coln
  DEGs <- degs_x
  remove(degs_x)
  # translate gene names in normalized data-file to human orthologs
  exprdata <- as.data.frame(exprdata)
  exprdata$'Mouse.gene.name' <- rownames(exprdata)
  exprdata <- inner_join(exprdata, orth)
  exprdata <- exprdata[!duplicated(exprdata$Gene.name),]
  rownames(exprdata) <- exprdata$Gene.name
  exprdata$Gene.name <- NULL
  exprdata$Mouse.gene.name <- NULL
  exprdata <- as.matrix(exprdata)
  remove(find_orthologs_mouse_to_human, list_orthologs, coln)
  
  # ANALYSIS
  #########################################################################################################
  
  # Analysis for all cell types
  
  Info<-list()
  counter=1
  
  cc <- unique(colnames(DEGs))
  all_exp_genes <- matrix(NA, nrow = nrow(exprdata), ncol = length(cc))
  # colnames(all_exp_genes) <- paste("Cluster_", cc, sep="")
  colnames(all_exp_genes) <- cc
  # s <- 1
  # t <- 24
  ntop_targets <- 250 # tested 700 as well, default == 250, last tested 600 did not work, 100o works
  for(s in 1:length(cc)){
    if(s == 1){
      rm(all_ligand_activity)
    }
    # sender cell type
    s_cells <- cell_types[cell_types[,colnames(cell_types)=='CellType_tissue']==cc[s],1]
    for(t in 1:length(cc)){
      print(paste(cc[s], " to ", cc[t], sep=""))
      # target cell type
      t_cells <- cell_types[cell_types[,colnames(cell_types)=='CellType_tissue'] == cc[t],1]
      
      # NICHENET ANALYSIS
      # define expressed genes
      if (length(s_cells) < 3 | length(t_cells) < 3){
        print('Less than 3 cells in one or both of the groups. Skip to next')
        next
      }
      unique(sort(paste(sapply(strsplit(colnames(exprdata), '_'), '[[', 1),
                        sapply(strsplit(colnames(exprdata), '_'), '[[', 2),
                        sapply(strsplit(colnames(exprdata), '_'), '[[', 3), sep = '_')))
      # which(colnames(exprdata) %in% s_cells)
      # which(s_cells %in% t_cells)
      exp_genes <- define_expressed_genes(exprdata, s_cells, t_cells, exp_cutoff = 1e-5)
      # define interesting genes
      int_genes <- genes_of_interest(sender_degs = DEGs[!is.na(DEGs[,colnames(DEGs) == cc[s]]),colnames(DEGs) == cc[s]],
                                     target_degs = DEGs[!is.na(DEGs[,colnames(DEGs) == cc[t]]),colnames(DEGs) == cc[t]],
                                     exp_genes)
      # identify ligand-receptor activity
      if(!all(is.na(int_genes[,1])) && !all(is.na(int_genes[,2]))){ # only analyse if cells express DEGs
        ligand_activity <- nichenet_analysis(exp_genes = exp_genes, int_genes = int_genes, ligand_target = ligand_target)
        ligand_activity2 <- ligand_activity[[2]]
        ligand_activity <- ligand_activity[[1]]
        
        Info[[counter]] <- ligand_activity2$info
        counter <- counter+1
        
        # Add gene target information
        ligand_activity$'target' <- NA
        ligand_activity$'target_weight' <- NA
        for (i in 1:length(ligand_activity$test_ligand)){
          ligand <- ligand_activity$test_ligand[i]
          active_ligand_target_links_df = ligand %>% lapply(get_weighted_ligand_target_links, geneset = int_genes[,2], ligand_target_matrix = ligand_target, n = ntop_targets) %>% bind_rows()
          ligand_activity$target[i] <- paste(active_ligand_target_links_df$target, collapse = '/')
          ligand_activity$target_weight[i] <- paste(active_ligand_target_links_df$weight, collapse = '/')
        }
        
        strsplit(as.character(ligand_activity[1,5]), '/')
        
        # write.table(ligand_activity$ranking, file = paste(dir.out, '/', cc[s], "_to_", cc[t], ".txt", sep=""),
        #             sep="\t", col.names = T, row.names = F)
        write.table(ligand_activity, file = paste(dir.out, '/', cc[s], "_to_", cc[t], ".txt", sep=""),
                    sep="\t", col.names = T, row.names = F)
        
        # comment below if predict_ligand_activities2
        if(!exists("all_ligand_activity")){
          all_ligand_activity <- cbind(ligand_activity, rep(cc[s], times = nrow(ligand_activity)), rep(cc[t], times = nrow(ligand_activity)))
          colnames(all_ligand_activity) <- c(colnames(ligand_activity), "Sender", "Target")
        } else {
          ligand_activity_x <- cbind(ligand_activity, rep(cc[s], times = nrow(ligand_activity)), rep(cc[t], times = nrow(ligand_activity)))
          colnames(ligand_activity_x) <- c(colnames(ligand_activity), "Sender", "Target")
          # all_ligand_activity <- rbind(all_ligand_activity, cbind(ligand_activity_x, rep(cc[s], times = nrow(ligand_activity_x)),
          #                                                         rep(cc[t], times = nrow(ligand_activity_x))))
          all_ligand_activity <- rbind(all_ligand_activity, ligand_activity_x)
        }
        # stop comment
        
      } else {
        print("No interesting DEGs")
      }
      
      # once for every cell type save expressed genes
      if(s == t){ # save expressed genes
        exp_genes <- exp_genes[!is.na(exp_genes[,1]),1]
        all_exp_genes[1:length(exp_genes),s] <- exp_genes
      }
    }
    
    all_ligand_activity <- all_ligand_activity %>% arrange(-pearson)
    write.table(all_ligand_activity, file = paste(dir.out, "/all_ligand_activity.txt", sep = ''), sep="\t", col.names = T, row.names = F)
    
  }
  
}




