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

scVI_v0.7.1.py

DEG_sort_significant.R

### Identification of cellular interactions and MCDM/MO-MCDM construction

To identify cell-cell interactions, we used NicheNet, a tool to model intercellular communications 
developed by [Browaeys et al., 2020](https://doi.org/10.1038/s41592-019-0667-5).
To define the genes of interest for identification of the intercellular interactions, we used the lists of DEGs for each cell type in each organ. 
The interactions were thus predicted, between each pair of cell types and each organ separetly. 

As the NicheNet database consists of intercellular interactions in humans, the gene names were translated to their human orthologs. 
The translation file used for our studies is provided in /data/orthologous_translation_file.txt. 
To run this code on any other dataset, it is recommended to download the orthologs from Ensembl as described [here](https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/).

The input files to run NicheNet_analysis.R, are;  
* the normalized expression data. Output from [scVI_v0.7.1.py](#data-normalization-and-differential-expression-analysis)
* the cell type identification file. Output from [cell type analysis script](#clustering-and-cell-type-identification)
* a translation file of human orthologs
* a matrix listing all DEGs for each cell type and organ combination, with celltype_organ over columns, 
as output from the [differential expression analysis](#data-normalization-and-differential-expression-analysis), DEG_sort_significant.R. 

The code output one tab separated txt file per interacting cell type pair, within and between organs, 
containing "test_ligand" (i.e upstream regulator of interaction), "auroc", "aupr", "pearson" (Pearson Correlation Coefficient, PCC), "target" (genes), and "target_weight" over the columns and interactions over rows. 

Additionally, it creates one file containing all predicted interactions between each pair of cell types, all_ligand_activity.txt, 
which apart from above mentioned columns contains "Sender" and "Target" cell type and organ (named: celltype_organ)

Note: that the interactions between organs in these files are not yet curated only to include URs secreated in blood

To curate the inter-organ interactions, only to include interactions through URs secreted in blood, 
and to sort out the strongest interactions (PCC > 0) considered for further analyses
we ran NicheNet_network_curation.R. As input to this script is the all_ligand_activity.txt output from NicheNet_analysis.R, 
and a curation file (eg., data/IPA/curation/curation_file.txt) containing information about the cellular location of each UR. 

The curation file can be retrieved from IPA by following [IPA - Generate curation file](#generate-curation-file)

It is important that the curation file contains two columns named "Symbol" (containing the gene names of the URs) and "Location" (containing information about cellular location). 
The URs which will be kept for inter-organ interactions by the code are those with "Location" == "Extracellular Space". 

The output from this script consists of one file containing all inter- and intra-organ interactions, all_pos_curated_ligand_activity.txt (Fig 1)

<img src="temp/images/output_nichenet_network_curation.PNG">

### Ranking of URs based on their downstream effect

/codes/rank_by_targets_and_heatmap.R

### Connective Pathway Analysis

## Meta analysis of 11 IMIDs

### Differential expression analysis

### ...

## Ingenuity pathway analysis

### Generate curation file

### Pathway analysis

### Prediction of upstream regulators (URs)
