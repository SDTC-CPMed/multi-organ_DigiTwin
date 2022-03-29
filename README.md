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

### Differential expression analysis

/codes/scVI_v0.7.1.py

### Identification of cellular interactions and MCDM/MO-MCDM construction

/codes/NicheNet_analysis.R

### Ranking of URs based on their downstream effect

/codes/rank_by_targets_and_heatmap.R

### Ingenuity Pathway Analysis for identification of enriched pathways

### Connective Pathway Analysis

## Meta analysis of 11 IMIDs

### Differential expression analysis

### ...
