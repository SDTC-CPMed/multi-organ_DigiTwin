# Project Title

## Instructions for running the codes

## scRNA-seq analysis for construction of MCDMs, MO-MCDMs, and Connective Pathway Analysis

The expected data inputs are five csv.gz tables (one per organ) with cells in columns and genes in rows. 

The data can be found at ...

### Quality assessment and full matrix construction

To ensure good quality data for downstream analyses, it is recommended to remove poor quality cells and genes from the analyses.
We performed these analyses for each organ separetly. 
We define and keep the good quality cells as those having a minimum of 400 transcripts, 200 genes, and less than 20% mitochondrial genes. 
The good quaity genes kept were defined as those being identified in at least 1% of the cells. 

Due to the risk of duplicates in the library resulting in two or more cells sharing a cell barcode, it is also recommend to remove outliers.
Based on empirical evaluation of the distributionan overestimation of transcripts count over the cells, we removed all cells with; 
Lung: >2000 transcripts; Spleen: >6000 transcripts; Muscle and Skin: >7000 transcripts; Joint: No outliers removed.

The expected input is for this part is a matrix with genes in rows and cells in columns. 
The quality of the data and removal of the outliers is done using sc_data_quality_sorting.R

### Clustering and cell type identification

### Data normalization and Differential expression analysis

/codes/scVI_v0.7.1.py

### Identification of cellular interactions and MCDM/MO-MCDM construction

To identify cell-cell interactions, we used NicheNet, a tool to model intercellular communications 
developed by [Browaeys et al., 2020](https://doi.org/10.1038/s41592-019-0667-5).

To define the genes of interest for identification of the intercellular interactions, we used the lists of DEGs for each cell type in each organ. 
The interactions were thus predicted, between each pair of cell types within each organ separetly. 
All interactions between all different cell types were combined into a Multicellular Disease Model (MCDM), for each organ separetly. 
Interactions were additionally identified between each cell type between different organs. These interactions were curated to only include those 
through ligands secreated into the blood (explain detailed criteria). A Multi-organ MCDM (MO-MCDM) was then constructed including all interactions 
between all the organs. 

As the NicheNet database consists of intercellular interactions in humans, the gene names were translated to their human orthologs. 
The translation file used for our studies is provided in /data/orthologous_translation_file.txt. 
To run this code on any other dataset, it is recommended to download the orthologs from Ensembl as described [here](https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/).

/codes/NicheNet_analysis.R

### Ranking of URs based on their downstream effect

/codes/rank_by_targets_and_heatmap.R

### Ingenuity Pathway Analysis for identification of enriched pathways

### Connective Pathway Analysis

## Meta analysis of 11 IMIDs

### Differential expression analysis

### ...
