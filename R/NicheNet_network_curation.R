rm(list=ls())
library(dplyr)


# define directories and file input
dir.data <- './data/NicheNet_analysis'
curation_file <- './data/IPA/curation/curation_file.txt'
dir.out <- './data/NicheNet_analysis_curated'
# dir.home <- getwd()
# dir.data <- paste(dir.home, '/results/NicheNet_analysis/cluster_ids_fromOleg_06_29_CellType/scVI_change', sep = '')
# curation_file <- paste(dir.home, '/results/IPA/cluster_ids_fromOleg_06_29_CellType/scVI_change/results/inter-tissue_interaction_sorting/aextra2_fix.txt', sep = '')
# dir.out <- paste(dir.home, '/results/postNicheNet_analysis/cluster_ids_fromOleg_06_29_CellType/scVI_change', sep = '')
if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print('dir.out was created')
}

# Load the data
all_ligand_activity <- read.table(paste(dir.data, "/all_ligand_activity.txt", sep = ''), sep = '\t', header = T)
cur <- read.table(curation_file)
# cur <- cur[,3:4]
# colnames(cur) <- c('Symbol', 'Location')
# write.table(cur, curation_file, sep = '\t', row.names = F)
# head(cur)

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



