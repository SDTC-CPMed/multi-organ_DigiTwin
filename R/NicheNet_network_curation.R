# Sandra Lilja
# R version 4.0.4
#'
#' Curate the inter-organ interactions for those throgu genes located in extracelllular space 
#' 
#' @param all_ligand_activity The combined MCDM output from NicheNet_analysis containing all inter- and intra-organ interactions
#' @param cur The curation file, containing information about cellular location. Must contain the columns "Symbol" 
#' (with gene symbols matching the Symbols in all_ligand_activity) and "Location" (where genes in "Extracellular Space" 
#' will be kept for potential inter-organ interactions)
#' @param dir.out  output directory
#' 
#' @export
#'

library(dplyr)

NicheNet_network_curation <- function(all_ligand_activity, cur, dir.out){
  # Ensure that all genes in interaction file are included in curation file
  if (length(which(unique(sort(all_ligand_activity$test_ligand)) %in% cur$Symbol == F))!= 0){
    print('Warning: Not all ligands are in curation file. Only the genes included in this file will be considered')
  } 
  # Keep the URs found in extracellular space
  cur <- cur[which(cur$Location %in% c('Extracellular Space')),]
  
  # include only interactions with positive PCC score
  all_ligand_activity <- all_ligand_activity[which(all_ligand_activity$pearson > 0),]
  
  # curate the inter-organ interactions only to contain URs from curation file
  # Prepare the data
  all_ligand_activity$'Sender_tissue' <- sapply(strsplit(as.character(all_ligand_activity$Sender), '_'), '[[', 2)
  all_ligand_activity$'Target_tissue' <- sapply(strsplit(as.character(all_ligand_activity$Target), '_'), '[[', 2)
  # Make one file including possible inter-tissue interactions only
  all_ligand_activity_inter <- all_ligand_activity[which(all_ligand_activity$Sender_tissue != all_ligand_activity$Target_tissue),]
  all_ligand_activity_inter <- all_ligand_activity_inter[which(all_ligand_activity_inter$test_ligand %in% cur$Symbol),]
  # combine with previous file
  all_ligand_activity_intra <- all_ligand_activity[which(all_ligand_activity$Sender_tissue == all_ligand_activity$Target_tissue),]
  all_ligand_activity <- full_join(all_ligand_activity_intra, all_ligand_activity_inter)
  
  # write to output
  write.table(all_ligand_activity, paste(dir.out, '/all_pos_curated_ligand_activity.txt', sep = ''), sep = '\t', col.names = T, row.names = F)
}



