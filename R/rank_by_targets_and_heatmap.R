# Sandra Lilja
# R version 4.0.4
#'
#' Rank the URs based on their downstream effect, and create heatmap of the ranked URs. 
#' 
#' @param 
#' @param 
#' @param 
#' 
#' @export
#'


library(dplyr)
library("pheatmap")

ints <- 'intra-tissue'
# ints <- 'inter-tissue'
# ints <- 'inter-and-intra-tissue'

dir.home <- getwd()
dir.target <- paste(dir.home, '/results/URs_prioritization/cluster_ids_fromOleg_06_29_CellType/scVI_change/centrality_analysis_all/data_in', sep = '')
dir.out <- paste(dir.home, '/results/rank_by_targets', sep = '')
if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print('dir.out was created')
}

# Load target information
target_list <- list.files(dir.target, pattern = 'ligand_activity')
length(target_list) == 1
targets <- read.table(paste(dir.target, target_list, sep = '/'), sep = '\t', header = T)
# Remove targets of interactions with negative PCC
targets <- targets[targets$pearson > 0,]

# subset, if only intra-tissiu interactions
if (ints == 'intra-tissue'){
  targets <- targets[which(sapply(strsplit(targets$Sender, '_'), '[[', 2) == 
            sapply(strsplit(targets$Target, '_'), '[[', 2)),]
} else if (ints == 'inter-tissue'){
  targets <- targets[which(sapply(strsplit(targets$Sender, '_'), '[[', 2) != 
                             sapply(strsplit(targets$Target, '_'), '[[', 2)),]
}

targets_backup <- targets

# define if only to rank on one tissue
Tis <- 'all'
# Tis <- 'Muscle'
# Tis <- 'Joint'
# Tis <- 'Lung'
# Tis <- 'Skin'
# Tis <- 'Spleen'

if (Tis != 'all'){
  targets <- targets_backup[grep(Tis, targets_backup$Sender),]
  targets <- targets[grep(Tis, targets$Target),]
}

# paste(sapply(strsplit(unique(sort(targets[targets$test_ligand == 'TNF','Target'])), '_'), '[[', 1), collapse = ', ')
# paste(sapply(strsplit(unique(sort(targets[targets$test_ligand == 'IL1B','Target'])), '_'), '[[', 1), collapse = ', ')
# paste(sapply(strsplit(unique(sort(targets[targets$test_ligand == 'TGFB1','Target'])), '_'), '[[', 1), collapse = ', ')

# List the URs
URs <- unique(sort(targets$test_ligand))

# Create N target matrix
N_target_data <- matrix(NA, ncol = length(URs), nrow = length(unique(sort(targets$Target)))) 
colnames(N_target_data) <- URs
rownames(N_target_data) <- unique(sort(targets$Target))
for (col in 1:length(colnames(N_target_data))){
  if (length(grep('\\.', col/10)) == 0){
    print(paste('column', as.character(col), 'of', as.character(length(colnames(N_target_data)))))
  }
  UR <- colnames(N_target_data)[col]
  for (row in 1:length(rownames(N_target_data))){
    CT <- rownames(N_target_data)[row]
    # sub_targets <- targets[grep(CT, targets$Target),]
    # sub_targets <- sub_targets[grep(UR, sub_targets$test_ligand),]
    sub_targets <- targets[which(targets$Target == CT),]
    sub_targets <- sub_targets[which(sub_targets$test_ligand == UR),]
    if (length(rownames(sub_targets)) == 0){
      N_target_data[row,col] <- 0
    } else {
      N_target_data[row,col] <- length(strsplit(sub_targets$target, '/')[[1]])
    }
  }
}

# rank order the URs based on: 
# 2. nr of unique target cell type and tissue combination, 
# 1. total nr of target genes, cell types, and tissue combinations
N_target_data <- N_target_data[,order(match(colnames(N_target_data), names(sort(colSums(N_target_data == 0)))))]
N_target_data <- N_target_data[,order(match(colnames(N_target_data), names(sort(colSums(N_target_data), decreasing = T))))]

# order the rows based on names
N_target_data <- N_target_data[order(sapply(strsplit(rownames(N_target_data), '_'), '[[', 1)),]
N_target_data <- N_target_data[order(sapply(strsplit(rownames(N_target_data), '_'), '[[', 2)),]

# create the heatmap for FC data
colSums(N_target_data)
range(colSums(N_target_data))
median(colSums(N_target_data))

bk2 <- c(0.001,seq(0.1,50,by=0.1))

my_palette <- c("white",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

# annotation_colors <- data.frame(row.names = rownames(N_target_data), 
#                                 Tissue = sapply(strsplit(rownames(N_target_data), '_'), '[[', 2),
#                                 CellType = sapply(strsplit(rownames(N_target_data), '_'), '[[', 1))
# 
# ann_colors = list(
#   Tissue = c(Joint = "#66C2A5", Lung = "#FC8D62", Muscle = "#8DA0CB", 
#              Skin = "#E78AC3", Spleen = "#FFD92F"),
#   CellType = c('B-cells' = "#A6CEE3", 'Dendritic-cells' = "#1F78B4", 
#                'Endothelial-cells' = "#B2DF8A", Erythrocytes = "#33A02C", 
#                Fibroblasts = "#FB9A99", Granulocytes = "#E31A1C",
#                Macrophages = "#FDBF6F", Monocytes = "#FF7F00", 
#                'NK-cells' = "#CAB2D6", 'T-cells' = "#6A3D9A", 
#                Unknown0 = "#FFFF99", Unknown16 = "#B15928",
#                Unknown6 = "#666666")
# )

ann_colors_in <- read.table(list.files(paste(dir.home, '/results', sep = ''), pattern = 'CellType_and_Tissue_colors', full.names = T), header = T)
annotation_colors <- data.frame(row.names = rownames(N_target_data), 
                                Tissue = sapply(strsplit(rownames(N_target_data), '_'), '[[', 2),
                                CellType = sapply(strsplit(rownames(N_target_data), '_'), '[[', 1))

# fix Unknown names
annotation_colors$'rn' <- rownames(annotation_colors)
annotation_colors <- left_join(annotation_colors, ann_colors_in, by = c('CellType' = 'Row'))
rownames(annotation_colors) <- annotation_colors$rn
annotation_colors$rn <- NULL
annotation_colors$CellType_IDnr <- NULL
annotation_colors$color <- NULL
annotation_colors$CellType <- NULL
colnames(annotation_colors) <- c('Tissue', 'CellType')

# ann_colors_in
Tissue <- ann_colors_in[1:5,2]
names(Tissue) <- ann_colors_in[1:5,1]
# Tissue <- Tissue[names(Tissue) %in% annotation_colors$Tissue]
CellType <- ann_colors_in[6:18,2]
names(CellType) <- ann_colors_in[6:18,3]
# CellType <- CellType[names(CellType) %in% annotation_colors$CellType]
ann_colors <- list(Tissue = Tissue, CellType = CellType)


ph <- pheatmap(N_target_data, color = my_palette, breaks = bk2, 
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         annotation_row = annotation_colors,
         annotation_colors = ann_colors)
ph

if (Tis == 'all'){
  write.csv(N_target_data, paste(dir.out, '/N_target-genes_matrix_', ints, '.csv', sep = ''))
  pdf(paste(dir.out, '/heatmap_N_target-genes_', ints, '.pdf', sep = ''), width = 14, height = 10)
  print(ph)
  dev.off()
} else {
  write.csv(N_target_data, paste(dir.out, '/N_target-genes_matrix_', ints, '_', Tis, '.csv', sep = ''))
  pdf(paste(dir.out, '/heatmap_N_target-genes_', ints, '_', Tis, '.pdf', sep = ''), width = 14, height = 10)
  print(ph)
  dev.off()
}

