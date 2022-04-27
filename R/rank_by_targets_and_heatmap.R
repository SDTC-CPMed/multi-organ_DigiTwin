# Sandra Lilja
# R version 4.0.4
#'
#' Rank the URs based on their downstream effect, and create heatmap of the ranked URs. 
#' 
#' @param targets The data frame containing all the interactions of interest, as output from NicheNet_network_curation
#' @param ints string 'inter' or 'intra' or 'all', dependent on which interactions to include for the ranking
#' @param dir.out output directory
#' @param Tis optional parameter, if ranking should be based on one tissue only, as named in targets$Target_tissue. Default is to include all. 
#' @param ann_colors_in optional parameter. Dataframe of preferred colors for labeling of tissues and cell types. Default = NULL
#' 
#' @export
#'


library(dplyr)
library("pheatmap")


rank_by_targets_and_heatmap <- function(targets, ints, dir.out, Tis = 'all', ann_colors_in = NULL){
  # subset, if only inter-/intra-tissue interactions
  if (ints == 'intra'){
    targets <- targets[which(targets$Sender_tissue == targets$Target_tissue),]
  } else if (ints == 'inter'){
    targets <- targets[which(targets$Sender_tissue != targets$Target_tissue),]
  }
  
  if (Tis != 'all'){
    targets <- targets[which(targets$Target_tissue == Tis),]
  }
  
  # List the URs
  URs <- unique(sort(targets$test_ligand))
  
  # Create N target matrix
  N_target_data <- matrix(NA, ncol = length(URs), nrow = length(unique(sort(targets$Target)))) 
  colnames(N_target_data) <- URs
  rownames(N_target_data) <- unique(sort(targets$Target))
  for (col in 1:length(colnames(N_target_data))){
    # if (length(grep('\\.', col/10)) == 0){
    #   print(paste('column', as.character(col), 'of', as.character(length(colnames(N_target_data)))))
    # }
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
  
  annotation_colors <- data.frame(row.names = rownames(N_target_data), 
                                  Tissue = sapply(strsplit(rownames(N_target_data), '_'), '[[', 2),
                                  CellType = sapply(strsplit(rownames(N_target_data), '_'), '[[', 1))
  # clean up cell type names
  annotation_colors$CellType <- gsub('-', ' ', annotation_colors$CellType)
  
  if (is.data.frame(ann_colors_in)){
    Tissue <- ann_colors_in[1:5,2]
    names(Tissue) <- ann_colors_in[1:5,1]
    Tissue <- Tissue[names(Tissue) %in% annotation_colors$Tissue]
    CellType <- ann_colors_in[6:18,2]
    names(CellType) <- ann_colors_in[6:18,3]
    CellType <- CellType[names(CellType) %in% annotation_colors$CellType]
    ann_colors <- list(Tissue = Tissue, CellType = CellType)
    ph <- pheatmap(N_target_data, color = my_palette, breaks = bk2, 
                   cluster_cols = F, cluster_rows = F, show_rownames = F,
                   annotation_row = annotation_colors,
                   annotation_colors = ann_colors)
  } else {
    ph <- pheatmap(N_target_data, color = my_palette, breaks = bk2, 
                   cluster_cols = F, cluster_rows = F, show_rownames = F,
                   annotation_row = annotation_colors)
  }
  
  if (Tis == 'all'){
    write.csv(N_target_data, paste(dir.out, '/UR-ranking_N_target-genes_matrix_', ints, '-interactions.csv', sep = ''))
    pdf(paste(dir.out, '/heatmap_UR-ranking_', ints, '-interactions.pdf', sep = ''), width = 14, height = 10)
    print(ph)
    dev.off()
  } else {
    write.csv(N_target_data, paste(dir.out, '/UR-ranking_N_target-genes_matrix_', ints, '-interactions_', Tis, '.csv', sep = ''))
    pdf(paste(dir.out, '/heatmap_UR-ranking_', ints, '-interactions_', Tis, '.pdf', sep = ''), width = 14, height = 10)
    print(ph)
    dev.off()
  }
  
}

