rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

# define directories
dir.home <- getwd()
dir.data <- paste(dir.home, '/results/allTissues_tissue-sample-BatchRemoval', sep = '')
dir.clust <- paste(dir.home, '/results/clusters_final_out', sep = '')
dir.out <- paste(dir.home, '/results/clusters_final_out', sep = '')

# read in files
print('load the data')
clusts <- read.csv(list.files(dir.clust, pattern = '^cluster_ids_fromOleg_', full.names = T))

exprdata <- read.csv(list.files(dir.data, pattern = 'normalized', full.names = T), row.names = 1)

# Define marker genes and cell types in specified order 
marker_genes <- c('Ms4a1', 'Ptprc', 'Cd22', 'Cd19',
                  'Cebpa', 'Nfe2', 'Ifit1', 'Cd7', 'Siglech',
                  'Cdh5', 'Pecam1', 'Fabp4',
                  'Alas2', 'Bpgm', 'Car2', 'Cpox', 'Hba.a1', 'Hbb.bt', 'Alad', 'Tfrc', 'Hba.a2',
                  'Col1a2', 'Col3a1', 'Fbln2', 'Fstl1', 'Gsn', 'Mmp2', 'Sparc', 'Vim',
                  'Camp', 'Ngp', 'S100a8', 'S100a9',
                  'Ly6g', 'Mmp9', 'Ltf',
                  'Cd14', 'C1qb', 'C1qc', 'Cd68', 'Ctss', 'Cxcl2', 'Gatm', 'Lyz2', 'Rgs1',
                  'F13a1', 'Ccl6', 'Cxcl3', 'Il1b', 'Ccr2',
                  'Klrk1', 'Il18rap', 'Gzma',
                  'Bcl11b', 'Il7r', 'Itk', 'Lef1',
                  'Bgn', 'Rgs16', 'Palld', 'Rgs5', 'Cpe', 'Sgce', 'Fhl1', 'Cd151', 'Adi1', 'Adgrf5',
                  'Krt14', 'Krt5', 'Mt2', 'BC100530', 'Fosl1', 'Ptgs2')

cellTypes <- unique(sort(clusts$CellType))
cellTypes <- cellTypes[c(1:6,12,7:11,13)]

# # Create seurat object
# CIA_seurat <- CreateSeuratObject(counts = t(exprdata), project = "CIA")
# CIA_seurat$Idents <- clusts$CellType
# 
# VlnPlot(CIA_seurat, marker_genes, idents = c('Joint', 'Spleen'), stack = TRUE, sort = TRUE) +
#   theme(legend.position = "none")


# create the data table for violin plots 
data_in <- clusts
data_in$CellTypeGuessed <- NULL
data_in_2 <- exprdata[,colnames(exprdata) %in% marker_genes]
data_in_2$cell_ID <- rownames(data_in_2)
data_in_2 <- reshape2::melt(data_in_2)
head(data_in)
head(data_in_2)
data_in <- full_join(data_in, data_in_2)

remove(data_in_2)

# # Define cell names and colors to use
ann_colors_in <- read.table(list.files(paste(dir.home, '/results', sep = ''), pattern = 'CellType_and_Tissue_colors', full.names = T), header = T)
ann_colors_in$Row <- gsub('-', ' ', ann_colors_in$Row)
ann_colors_in <- ann_colors_in[ann_colors_in$Row %in% cellTypes,]
ann_colors_in <- ann_colors_in[match(cellTypes, ann_colors_in$Row),]

# data_in <- left_join(data_in, ann_colors_in, by = c('CellType' = 'Row'))
# head(data_in)

# Update Feature and Identity factor orders
data_in$CellType <- factor(data_in$CellType, levels = cellTypes)
data_in$variable <- factor(data_in$variable, levels = marker_genes)
data_in$CellTypes_fix

# add and order corrected cell type names
data_in <- left_join(data_in, ann_colors_in, by = c('CellType' = 'Row'))
data_in$CellTypes_fix <- factor(data_in$CellTypes_fix, levels = ann_colors_in$CellTypes_fix)

head(data_in)

b <- ggplot(data_in, aes(value, factor(CellTypes_fix), fill = variable)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE,
              aes(fill = factor(CellTypes_fix))) +
  scale_fill_manual(values = ann_colors_in$color) +
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


b

# write to output
outfile <- paste('Violin_markerGenes_vs_CellTypes.pdf', sep = '')
pdf(paste(dir.out, '/', outfile, sep = ''), width = 22, height = 5)
b
dev.off()


################
## boxplot
## % cell types over samples
#### Plot the distribution of different cell types for each separate timepoint and treatment
detach('package:dplyr')
library(plyr)
library(scales)
library(reshape2)
library(dplyr)

head(clusts)
clusts$'Mouse_ID' <- sapply(strsplit(clusts$cell_ID, '_'), '[[', 3)

### Create data frame 'Mdf' with unique rows containing column 'value'
### specifying the ratio of the cell type in its certain timepoint and condition.
df <- ddply(clusts,.(CellType, disease_group, tissue, Mouse_ID),nrow)
colnames(df) <- c('CellType', 'disease_group', 'tissue', 'Mouse_ID', 'ratio')
df <- left_join(df, ann_colors_in, by = c('CellType' = 'Row'))
df$disease_group <- gsub('Sick', 'CIA', df$disease_group)

# plot % cells of each cell type in each tissue and disease stage
newplot <- ggplot(df, aes(x=disease_group, y=ratio, fill=CellTypes_fix)) +
  geom_bar(stat="identity", position = "fill",
           aes(fill = CellTypes_fix), 
           width = 0.95) +
  scale_fill_manual(values = ann_colors_in$color[order(ann_colors_in$CellTypes_fix)]) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(~tissue) + theme_bw() +
  theme(
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(2, 'mm')
  )

newplot

outfile <- paste('Proportions_CellTypes_per_tissue_and_state.pdf', sep = '')
pdf(paste(dir.out, '/', outfile, sep = ''))
newplot
dev.off()


# plot % cells of each disease stage in each cell type
newplot <- ggplot(df, aes(x=ratio, y=CellTypes_fix, fill=disease_group)) +
  geom_bar(stat="identity", position = "fill",
           aes(fill = disease_group), 
           width = 0.95) +
  scale_fill_manual(values = c('indianred1', 'lightblue3')) +
  scale_x_continuous(labels = percent_format()) +
  theme(
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(2, 'mm'),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA)
  )

newplot

outfile <- paste('Proportions_state_per_CellType.pdf', sep = '')
pdf(paste(dir.out, '/', outfile, sep = ''))
newplot
dev.off()


# plot % cells of each tissue in each cell type
ann_colors_in_2 <- read.table(list.files(paste(dir.home, '/results', sep = ''), pattern = 'CellType_and_Tissue_colors', full.names = T), header = T)
ann_colors_in_2$Row <- gsub('-', ' ', ann_colors_in_2$Row)
tissues <- unique(sort(df$tissue))
ann_colors_in_2 <- ann_colors_in_2[ann_colors_in_2$Row %in% tissues,]
ann_colors_in_2$CellType_IDnr <- NULL

newplot <- ggplot(df, aes(x=ratio, y=CellTypes_fix, fill=tissue)) +
  geom_bar(stat="identity", position = "fill",
           aes(fill = tissue), 
           width = 0.95) +
  scale_fill_manual(values = ann_colors_in_2$color) +
  scale_x_continuous(labels = percent_format()) +
  theme(
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(2, 'mm'),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA)
  )

newplot

outfile <- paste('Proportions_tissues_per_CelType.pdf', sep = '')
pdf(paste(dir.out, '/', outfile, sep = ''))
newplot
dev.off()

# plot % cells from each tissue in each sample
newplot <- ggplot(df, aes(x=Mouse_ID, y=ratio, fill=tissue)) +
  geom_bar(stat="identity", position = "fill",
           aes(fill = tissue), 
           width = 0.95) +
  scale_fill_manual(values = ann_colors_in_2$color) +
  scale_y_continuous(labels = percent_format()) +
  theme(
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(2, 'mm')
  )

newplot

outfile <- paste('Proportions_tissue_per_MouseID.pdf', sep = '')
pdf(paste(dir.out, '/', outfile, sep = ''))
newplot
dev.off()

# plot % cells of each cell type in each sample
newplot <- ggplot(df, aes(x=Mouse_ID, y=ratio, fill=CellTypes_fix)) +
  geom_bar(stat="identity", position = "fill",
           aes(fill = CellTypes_fix), 
           width = 0.95) +
  scale_fill_manual(values = ann_colors_in$color[order(ann_colors_in$CellTypes_fix)]) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(~tissue) + theme_bw() +
  theme(
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(2, 'mm')
  )

newplot

outfile <- paste('Proportions_CellTypes_per_MouseID_and_tissue.pdf', sep = '')
pdf(paste(dir.out, '/', outfile, sep = ''), width = 14)
newplot
dev.off()




