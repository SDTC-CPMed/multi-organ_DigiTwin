# Sandra Lilja
# R version 4.0.4
#'
#' Plot the expression level and/or FC of a set of genes 
#' 
#' @param exprdata The normalized single-cell expression matrix
#' @param DE_data The lists of tables including significant DEGs
#' @param clusts The clusters and cell type information, containing the column 'CellType'
#' @param marker_genes  a vector containing the marker genes of interest
#' 
#' @export
#'

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(pheatmap)

marker_gene_expression <- function(exprdata, clusts, marker_genes){
  # gsub to match gene symbols
  colnames(exprdata) <- gsub('\\.', '-', colnames(exprdata))
  
  cellTypes <- unique(sort(clusts$CellType))
  # cellTypes <- cellTypes[c(1:6,12,7:11,13)]
  
  # create the data table for violin plots 
  data_in <- clusts
  data_in_2 <- exprdata[,colnames(exprdata) %in% marker_genes]
  data_in_2$cell_ID <- rownames(data_in_2)
  data_in_2 <- reshape2::melt(data_in_2)
  # head(data_in)
  # head(data_in_2)
  data_in <- full_join(data_in, data_in_2)
  remove(data_in_2)
  
  # # Define cell names and colors to use
  # ann_colors_in <- read.table(list.files(paste(dir.home, '/results', sep = ''), pattern = 'CellType_and_Tissue_colors', full.names = T), header = T)
  # ann_colors_in$Row <- gsub('-', ' ', ann_colors_in$Row)
  # ann_colors_in <- ann_colors_in[ann_colors_in$Row %in% cellTypes,]
  # ann_colors_in <- ann_colors_in[match(cellTypes, ann_colors_in$Row),]
  
  # Update Feature and Identity factor orders
  data_in$CellType <- factor(data_in$CellType, levels = cellTypes)
  data_in$variable <- factor(data_in$variable, levels = marker_genes)
  data_in$CellTypes_fix
  
  # # add and order corrected cell type names
  # data_in <- left_join(data_in, ann_colors_in, by = c('CellType' = 'Row'))
  # data_in$CellTypes_fix <- factor(data_in$CellTypes_fix, levels = ann_colors_in$CellTypes_fix)
  
  head(data_in)
  
  b <- ggplot(data_in, aes(value, factor(CellType), fill = variable)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE,
                aes(fill = factor(CellType))) +
    # scale_fill_manual(values = ann_colors_in$color) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(cols = vars(variable), scales = "free")  +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          panel.background = element_rect(fill = NA, color = "grey"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.title.y.left = element_blank())
  
  # # write to output
  # outfile <- paste('Violin_markerGenes_vs_CellTypes.pdf', sep = '')
  # pdf(paste(dir.clust, '/', outfile, sep = ''), width = 22, height = 5)
  # b
  # dev.off()
  
  return(b)
}

marker_gene_DE <- function(DE_data, dir.DE, clusts, marker_genes){
  # read in the DEGs
  DE_all <- c()
  for (i in 1:length(DE_data)){
    DE_all[[i]] <- read.csv(paste(dir.DE, DE_data[i], sep = '/'))
  }
  # head(clusts)

  # create matrix for heatmap
  cellTypes <- sapply(strsplit(DE_data, '_'), '[[', 3)
  names(DE_all) <- cellTypes
  M <- matrix(data = NA, nrow = length(cellTypes), ncol = length(marker_genes))
  colnames(M) <- marker_genes
  rownames(M) <- cellTypes
  
  for (i in 1:length(rownames(M))){
    cellType <- rownames(M)[i]
    DE_x <- DE_all[[which(names(DE_all) == cellType)]]
    for (j in 1:length(colnames(M))){
      gene <- colnames(M)[j]
      if (length(which(DE_x$X == gene)) == 0){
        M[i,j] <- 0
      } else if (length(which(DE_x$X == gene)) == 1){
        if (DE_x[which(DE_x$X == gene),'lfc_mean']>0){
          M[i,j] <- 1
        } else if (DE_x[which(DE_x$X == gene),'lfc_mean']<0){
          M[i,j] <- 0
        } else if (DE_x[which(DE_x$X == gene),'lfc_mean'] == 0){
          print(paste('Error:', i, j, 'gene DE but lfc == 0'))
        }
      } else if (length(which(DE_x$X == gene)) > 1)
        print(paste('Error:', i, j, 'gene identified', length(which(DE_x$X == gene)), 'times'))
    }
  }
  ph <- pheatmap(M, 
           cluster_rows = F,
           cluster_cols = F)
  return(ph)
}


# ################
# ## boxplot
# ## % cell types over samples
# #### Plot the distribution of different cell types for each separate timepoint and treatment
# detach('package:dplyr')
# library(plyr)
# library(scales)
# library(reshape2)
# library(dplyr)
# 
# head(clusts)
# clusts$'Mouse_ID' <- sapply(strsplit(clusts$cell_ID, '_'), '[[', 3)
# 
# ### Create data frame 'Mdf' with unique rows containing column 'value'
# ### specifying the ratio of the cell type in its certain timepoint and condition.
# df <- ddply(clusts,.(CellType, disease_group, tissue, Mouse_ID),nrow)
# colnames(df) <- c('CellType', 'disease_group', 'tissue', 'Mouse_ID', 'ratio')
# df <- left_join(df, ann_colors_in, by = c('CellType' = 'Row'))
# df$disease_group <- gsub('Sick', 'CIA', df$disease_group)
# 
# # plot % cells of each cell type in each tissue and disease stage
# newplot <- ggplot(df, aes(x=disease_group, y=ratio, fill=CellTypes_fix)) +
#   geom_bar(stat="identity", position = "fill",
#            aes(fill = CellTypes_fix), 
#            width = 0.95) +
#   scale_fill_manual(values = ann_colors_in$color[order(ann_colors_in$CellTypes_fix)]) +
#   scale_y_continuous(labels = percent_format()) +
#   facet_grid(~tissue) + theme_bw() +
#   theme(
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_blank(),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.spacing = unit(2, 'mm')
#   )
# 
# newplot
# 
# outfile <- paste('Proportions_CellTypes_per_tissue_and_state.pdf', sep = '')
# pdf(paste(dir.clust, '/', outfile, sep = ''))
# newplot
# dev.off()
# 
# 
# # plot % cells of each disease stage in each cell type
# newplot <- ggplot(df, aes(x=ratio, y=CellTypes_fix, fill=disease_group)) +
#   geom_bar(stat="identity", position = "fill",
#            aes(fill = disease_group), 
#            width = 0.95) +
#   scale_fill_manual(values = c('indianred1', 'lightblue3')) +
#   scale_x_continuous(labels = percent_format()) +
#   theme(
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_blank(),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.spacing = unit(2, 'mm'),
#     panel.background = element_blank(),
#     panel.border = element_rect(fill = NA)
#   )
# 
# newplot
# 
# outfile <- paste('Proportions_state_per_CellType.pdf', sep = '')
# pdf(paste(dir.clust, '/', outfile, sep = ''))
# newplot
# dev.off()
# 
# 
# # plot % cells of each tissue in each cell type
# ann_colors_in_2 <- read.table(list.files(paste(dir.home, '/results', sep = ''), pattern = 'CellType_and_Tissue_colors', full.names = T), header = T)
# ann_colors_in_2$Row <- gsub('-', ' ', ann_colors_in_2$Row)
# tissues <- unique(sort(df$tissue))
# ann_colors_in_2 <- ann_colors_in_2[ann_colors_in_2$Row %in% tissues,]
# ann_colors_in_2$CellType_IDnr <- NULL
# 
# newplot <- ggplot(df, aes(x=ratio, y=CellTypes_fix, fill=tissue)) +
#   geom_bar(stat="identity", position = "fill",
#            aes(fill = tissue), 
#            width = 0.95) +
#   scale_fill_manual(values = ann_colors_in_2$color) +
#   scale_x_continuous(labels = percent_format()) +
#   theme(
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_blank(),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.spacing = unit(2, 'mm'),
#     panel.background = element_blank(),
#     panel.border = element_rect(fill = NA)
#   )
# 
# newplot
# 
# outfile <- paste('Proportions_tissues_per_CelType.pdf', sep = '')
# pdf(paste(dir.clust, '/', outfile, sep = ''))
# newplot
# dev.off()
# 
# # plot % cells from each tissue in each sample
# newplot <- ggplot(df, aes(x=Mouse_ID, y=ratio, fill=tissue)) +
#   geom_bar(stat="identity", position = "fill",
#            aes(fill = tissue), 
#            width = 0.95) +
#   scale_fill_manual(values = ann_colors_in_2$color) +
#   scale_y_continuous(labels = percent_format()) +
#   theme(
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_blank(),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.spacing = unit(2, 'mm')
#   )
# 
# newplot
# 
# outfile <- paste('Proportions_tissue_per_MouseID.pdf', sep = '')
# pdf(paste(dir.clust, '/', outfile, sep = ''))
# newplot
# dev.off()
# 
# # plot % cells of each cell type in each sample
# newplot <- ggplot(df, aes(x=Mouse_ID, y=ratio, fill=CellTypes_fix)) +
#   geom_bar(stat="identity", position = "fill",
#            aes(fill = CellTypes_fix), 
#            width = 0.95) +
#   scale_fill_manual(values = ann_colors_in$color[order(ann_colors_in$CellTypes_fix)]) +
#   scale_y_continuous(labels = percent_format()) +
#   facet_grid(~tissue) + theme_bw() +
#   theme(
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_blank(),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.spacing = unit(2, 'mm')
#   )
# 
# newplot
# 
# outfile <- paste('Proportions_CellTypes_per_MouseID_and_tissue.pdf', sep = '')
# pdf(paste(dir.clust, '/', outfile, sep = ''), width = 14)
# newplot
# dev.off()




