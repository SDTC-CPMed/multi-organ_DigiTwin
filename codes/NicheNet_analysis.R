rm(list=ls())
library(nichenetr)
library(dplyr)

# samp <- 'cluster_ids_fromOleg'
samp <- 'cluster_ids_fromOleg_06_29_CellType'
# samp <- 'cluster_ids_fromOleg_06_29_CellTypeGuessed'
species <- 'mouse'
# ct <- 'cluster_ids_fromOleg.csv'
ct <- 'cluster_ids_fromOleg_06_29.csv'
# method = 'vanilla'
method = 'change'


# define directories
dir.home <- getwd()
# dir.data <- paste(dir.home, '/results/allTissues_noBatchRemoval/scVI_out_0.5', sep = '')
# dir.data <- paste(dir.home, '/results/allTissues_noBatchRemoval/scVI_out_1.5', sep = '')
# dir.data <- paste(dir.home, '/results/allTissues_noBatchRemoval', sep = '')
dir.data <- paste(dir.home, '/results/allTissues_tissue-sample-BatchRemoval', sep = '')
dir.ct <- paste(dir.home, '/results/clusters_final_out', sep = '')
dir.degs <- paste(dir.home, '/results/DEG_analysis/', samp, '/', method, '_mode/fdr_sorted', sep = '')
# dir.out <- paste(dir.home, '/results/NicheNet_analysis/', samp, sep = '')
# dir.orth <- paste(dir.home, '/results/orthologous_translation/', samp, sep = '')
dir.out <- paste(dir.home, '/results/NicheNet_analysis/', samp, '/scVI_change', sep = '')
# dir.out <- paste(dir.home, '/results/NicheNet_analysis/', samp, '/scVI_change/targets_fix_700', sep = '')
dir.orth <- paste(dir.home, '/results/orthologous_translation/', samp, '/scVI_change', sep = '')
# dir.out <- paste(dir.home, '/results/NicheNet_analysis/', samp, '/scVI_change_test', sep = '')
# dir.orth <- paste(dir.home, '/results/orthologous_translation/', samp, '/scVI_change_test', sep = '')
if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print('dir.out was created')
}
if (dir.exists(dir.orth)==FALSE){
  dir.create(dir.orth, recursive = T)
  print('dir.orth was created')
}


# scVI adjusted scRNA-seq mouse data 
data <- read.csv(paste(dir.data, "/normalized_expression_matrix.csv", sep = ''), row.names = 1)
if (species == 'human'){
  data <- as.matrix(data)
}
data <- t(data)
data[1:5,1:5]
# mode(data) <- "numeric"

# Cell types
cell_types <- as.matrix(read.csv(paste(dir.ct, ct, sep = '/')))
# cell_types[,2] <- as.numeric(cell_types[,2])
# remove any unused column. CellType_tissue must be in the fifth column
head(cell_types)
if (samp == 'cluster_ids_fromOleg_06_29_CellType'){
  cell_types <- cell_types[,c(1:2,4:5)]
} else if (samp == 'cluster_ids_fromOleg_06_29_CellTypeGuessed'){
  cell_types <- cell_types[,c(1,3:5)]
  colnames(cell_types)[2] <- 'CellType'
}
# add CellType_tissue to new column
CellType_tissue <- as.vector(sub(' ', '_', paste(sub(' ', '-', cell_types[,'CellType']), cell_types[,'tissue'])))
cell_types <- cbind(cell_types, CellType_tissue)
remove(CellType_tissue)


# Get DEGs for cell types
degs <- as.matrix(read.csv(paste(dir.degs, "/DEGs_Sick_vs_Healthy.csv", sep = '')))
colnames(degs) <- gsub('\\.', '-', colnames(degs))

# Get ligand-target matrix
ligand_target <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds")) # targets = rows, ligands = columns

# Get ligand-receptor interactions
lr_network = as.matrix(readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds")))

# test
# signaling = as.matrix(readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds")))
# gr_network = as.matrix(readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds")))

# CONVERT GENES TO HUMAN ORTHOLOGS
#########################################################################################################
if (species == 'mouse'){
  list_orthologs  <- function(mouse_gene_list, orth){
    orth_sub = orth[which(orth$Mouse.gene.name %in% mouse_gene_list == TRUE),]
    # human_gene_orthologs <- as.character(orth$Gene.name[which(orth$Mouse.gene.name %in% mouse_gene_list == TRUE)])
    return(orth_sub)
  }
  find_orthologs_mouse_to_human  <- function(mouse_gene_list, orth){
    human_gene_orthologs <- as.character(orth$Gene.name[which(orth$Mouse.gene.name %in% mouse_gene_list == TRUE)])
    return(human_gene_orthologs)
  }
  # orth <- read.csv('/home/sanli71/fillager_OmikaHome/warefolder/data/Ensembl/GRCh38.p13.human_mouse_orthologs_ensembl_mart_export_200611.csv')
  orth <- read.csv(paste(dir.home, '/Ensembl/GRCh38.p13.human_mouse_orthologs_ensembl_mart_export_200611.csv', sep = ''))
  orth <- orth[,c('Gene.name', 'Mouse.gene.name')]
  orth <- distinct(orth)
  head(orth)
  for (col in 1:length(colnames(degs))){
    # col <- 8
    mouse_gene_list <- degs[,col]
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
    degs[,col] <- NA
    if (length(human_gene_orthologs) != 0){
      degs[1:length(human_gene_orthologs),col] <- human_gene_orthologs
    }
  }
  # remove samples where no orthologous genes could be found
  coln <- c()
  for (i in 1:length(colnames(degs))){
    if (length(unique(degs[,i])) == 1){
      if (is.na(unique(degs[,i]))){
        next
      }
    }
    else {
      coln <- c(coln, colnames(degs)[i])
      if (exists('degs_x')){
        degs_x <- cbind(degs_x, degs[,i])
      }
      else {
        degs_x <- degs[,i]
      }
    }
  }
  colnames(degs_x) <- coln
  degs <- degs_x
  remove(degs_x)
  # 
  orth_sub <- distinct(orth_sub)
  # data_backup <- data
  # data <- data_backup
  data <- as.data.frame(data)
  data$'Mouse.gene.name' <- rownames(data)
  data <- inner_join(data, orth)
  data <- data[!duplicated(data$Gene.name),]
  rownames(data) <- data$Gene.name
  data$Gene.name <- NULL
  data$Mouse.gene.name <- NULL
  data <- as.matrix(data)
  remove(orth, find_orthologs_mouse_to_human, coln)
}
write.table(orth_sub, paste(dir.orth, '/orthologous_translation_NicheNet_analysis.txt', sep = ''), sep = '\t', row.names = F)

# ANALYSIS
#########################################################################################################

# STEPS referred to below are from this vignette:
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md


# STEP 1: define expressed genes in sender and receiver cell population
# s = cell IDs sender
# t = cell IDs target
define_expressed_genes <- function(data, s = colnames(data), t = colnames(data), exp_cutoff = 0.1){
  
  # exp_cutoff <- 0.0001
  
  # reverse log10 tranformation
  ### Why ?
  expressed_genes_sender <- 10**data[,colnames(data) %in% s]-1
  expressed_genes_target <- 10**data[,colnames(data) %in% t]-1
  # expressed_genes_sender <- data[,colnames(data) %in% s]
  # expressed_genes_target <- data[,colnames(data) %in% t]
  # apply NicheNet recommended expression threshold
  temp <- log2(rowMeans(expressed_genes_sender)+1) >= exp_cutoff
  expressed_genes_sender <- rownames(expressed_genes_sender)[temp]
  temp <- log2(rowMeans(expressed_genes_target)+1) >= exp_cutoff
  expressed_genes_target <- rownames(expressed_genes_target)[temp]
  print(paste("n genes expressed in sender:", length(expressed_genes_sender), " n genes expressed in target:", length(expressed_genes_target), sep=" "))
  
  exp_genes <- cbind(c(expressed_genes_sender, rep(NA, times = nrow(data)-length(expressed_genes_sender))),
                     c(expressed_genes_target, rep(NA, times = nrow(data)-length(expressed_genes_target))))
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
  
  int_genes <- cbind(c(sender_interest_gene, rep(NA, times = nrow(degs)-length(sender_interest_gene))),
                     c(target_interest_gene, rep(NA, times = nrow(degs)-length(target_interest_gene))))
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




# Analysis for all cell types
# This calls all previously defined functions.

Info<-list()
counter=1

# c <- unique(cell_types[,5])

c <- unique(colnames(degs))
all_exp_genes <- matrix(NA, nrow = nrow(data), ncol = length(c))
# colnames(all_exp_genes) <- paste("Cluster_", c, sep="")
colnames(all_exp_genes) <- c
# s <- 10
# t <- 24
ntop_targets <- 250 # tested 700 as well, default == 250, last tested 600 did not work, 100o works
for(s in 1:length(c)){
  if(s == 1){
    rm(all_ligand_activity)
  }
  # sender cell type
  s_cells <- cell_types[cell_types[,5]==c[s],1]
  for(t in 1:length(c)){
    print(paste(c[s], " to ", c[t], sep=""))
    # target cell type
    t_cells <- cell_types[cell_types[,5] == c[t],1]
    
    # NICHENET ANALYSIS
    # define expressed genes
    if (length(s_cells) < 3 | length(t_cells) < 3){
      print('Less than 3 cells in one or both of the groups. Skip to next')
      next
    }
    unique(sort(paste(sapply(strsplit(colnames(data), '_'), '[[', 1),
                      sapply(strsplit(colnames(data), '_'), '[[', 2),
                      sapply(strsplit(colnames(data), '_'), '[[', 3), sep = '_')))
    which(colnames(data)[1:15] %in% s_cells)
    which(s_cells %in% t_cells)
    exp_genes <- define_expressed_genes(data, s_cells, t_cells, exp_cutoff = 1e-5)
    # define interesting genes
    int_genes <- genes_of_interest(sender_degs = degs[!is.na(degs[,colnames(degs) == c[s]]),colnames(degs) == c[s]],
                                   target_degs = degs[!is.na(degs[,colnames(degs) == c[t]]),colnames(degs) == c[t]],
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
      
      # write.table(ligand_activity$ranking, file = paste(dir.out, '/', c[s], "_to_", c[t], ".txt", sep=""),
      #             sep="\t", col.names = T, row.names = F)
      write.table(ligand_activity, file = paste(dir.out, '/', c[s], "_to_", c[t], ".txt", sep=""),
                  sep="\t", col.names = T, row.names = F)

      # comment below if predict_ligand_activities2
      if(!exists("all_ligand_activity")){
        all_ligand_activity <- cbind(ligand_activity, rep(c[s], times = nrow(ligand_activity)), rep(c[t], times = nrow(ligand_activity)))
        colnames(all_ligand_activity) <- c(colnames(ligand_activity), "Sender", "Target")
      } else {
        ligand_activity_x <- cbind(ligand_activity, rep(c[s], times = nrow(ligand_activity)), rep(c[t], times = nrow(ligand_activity)))
        colnames(ligand_activity_x) <- c(colnames(ligand_activity), "Sender", "Target")
        # all_ligand_activity <- rbind(all_ligand_activity, cbind(ligand_activity_x, rep(c[s], times = nrow(ligand_activity_x)),
        #                                                         rep(c[t], times = nrow(ligand_activity_x))))
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


#now computing correlations per liggand over all cell types

# # correlation ranking. 
# # NOTE: int_genes created per loop in above script. This only calcluates for latest interaction pair now
# PotentialLigands<-unique(int_genes[!is.na(int_genes[,1]),1])
# nLig=length(PotentialLigands)
# Correlations=numeric(nLig)
# 
# for (i in 1:nLig){
#   liggand=PotentialLigands[[i]]
#   rel=lapply(Info, function (x) x[[liggand]]$info)%>%bind_rows()
#   Correlations[i]=cor(rel$response, rel$prediction)
# }
# 
# #This is the table containing overall correlation rankings per ligand.
# overall_corr_table=data.frame(ligand=PotentialLigands, corr=Correlations)

# show histogram of ligand activity scores
library(ggplot2)
# ligand_activities <- as.data.frame(ligand_activity)
ligand_activities <- as.data.frame(all_ligand_activity)
ligand_activities$pearson <- as.numeric(as.character(ligand_activities$pearson))
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities  %>% top_n(5000, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities  %>% top_n(5, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

# Identify the best upstream ligands
best_upstream_ligands = ligand_activities %>% top_n(5000, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands = ligand_activities %>% top_n(5, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands_unique = as.character(unique(best_upstream_ligands))
head(best_upstream_ligands_unique)
length(best_upstream_ligands_unique)

# Infer target genes of top ranked ligands
# NOTE:
# int_genes not correct as it only contains the genes of interest for one of the interactions
# Redo: one analyses per interaction pair. Use best_upstream_ligands_unique which are present in interaction X and genes of interest for interaction X
# Save the genes of interest for each interaction from the above analysis
active_ligand_target_links_df = best_upstream_ligands_unique %>% lapply(get_weighted_ligand_target_links, geneset = int_genes, ligand_target_matrix = ligand_target, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
## [1] 1351
head(active_ligand_target_links_df)

# prepare for visualization
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target, cutoff = 0.25)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

# Visualize the ligand target activity
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network


'
# Overview - how many cell-cell interactions between cell X and cell Y
for(i in 1:length(unique(all_ligand_activity[,5]))){
if(i == 1){
rm(overview_ligand_activity)
}
for(j in 1:length(unique(all_ligand_activity[,6]))){
t_i <- unique(all_ligand_activity[,5])[i]
t_j <- unique(all_ligand_activity[,6])[j]

temp <- all_ligand_activity

if(sum(temp[,5]== t_i)>1){
temp <- temp[temp[,5]== t_i,]
if(sum(temp[,6]== t_j)>1){
temp <- temp[temp[,6]== t_j,]
temp <- nrow(temp)
} else {
if(sum(temp[,6]== t_j)==1){
temp <- 1
} else {
temp <- 0
}
}
} else {
if(sum(temp[,5]== t_i)==1){
if(which(temp[,5]==t_i) %in% which(temp[,6] == t_j)){
temp <- 1
} else {
temp <- 0
}
} else {
temp <- 0
}
}

if(!exists("overview_ligand_activity")){
overview_ligand_activity <- matrix(c(t_i, t_j, temp), nrow = 1)
colnames(overview_ligand_activity) <- c("Sender", "Target", "n_interactions")
} else {
overview_ligand_activity <- rbind(overview_ligand_activity, c(t_i, t_j, temp))
}
}
}
rm(t_i, t_j, i, j, temp)
mode(overview_ligand_activity) <- "numeric"
write.table(overview_ligand_activity, file = "Output/NicheNet/overview_interactions.txt", sep="\t", col.names = T, row.names = F)
'


# Used for choosing an arbitary threshold for when a gene is considered to be expressed in a given cell population
# out <- foreach(i = seq(from = 0, to = 0.2, by = 0.01), .combine = "cbind") %do% {
#   t_x <- vector()
#   for(j in 0:max(as.numeric(cell_types[,2]))){
#   # select cell population
#   temp <- data[,colnames(data)%in%cell_types[as.numeric(cell_types[,2])==j,1]]
#   # Criteria
#   if(ncol(temp)>0){
#   x <- rowMeans(temp>=i)
#   # Print result
#   t_x[j] <- length(x[x >= 0.1])
#   } else {
#   t_x[j] <- 0
#   }
#   }
#   return(t_x)
# }
# rownames(out) <- paste("Cluster_", 0:(nrow(out)-1), sep="")
# colnames(out) <- paste("Cutoff_", seq(from = 0, to = 0.2, by = 0.01), sep="")
# out <- out[rownames(out)%in% colnames(degs),]

