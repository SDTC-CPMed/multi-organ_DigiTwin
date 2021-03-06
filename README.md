# Multi-organ digital twin

## Analysis protocols

In this project, we started with multi-organ single cell analysis of a mouse model of collagen-induced arthritis. We analyzed five different organ samples (joint, lung, muscle, skin, and spleen) from 6 CIA and 4 healthy control mice. The expression matrix containing this data can be found in [GSE206651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206651) (GSE206651_UMI_expression_matrix.txt.gz, *will be made available upon publication*). We further analyzed the data following the steps in [CIA-analyses.md](./vignettes/CIA-analyses.md).

Next, we performed a meta analysis of 32 different datasets from variable locations/disease states of ten different immune-mediated inflammatory diseases (IMIDs), namely rheumatoid arthritis (RA), ulcerative colitis (UC), Crohn’s disease (CD), psoriasis (PSO), Sjögren’s syndrome (SS), systemic sclerosis (SSc), atopic dermatitis (AD), juvenile myositis (JM), ‘at risk for’ type 1 diabetes (T1D), and systemic lupus erythematosus (SLE), sampled from GSE16161, GSE32924, GSE16879, GSE179285, GSE75214, GSE81071, GSE148810, GSE112943, GSE32591, GSE14905, GSE181318, GSE1919, GSE55235, GSE176510, GSE40568, GSE81292, GSE95065, GSE66413, and GSE11223. For each of the datasets, the DEGs were calculated using GEO2R, between sick and healthy individuals, as further described by Lilja et al., *in manuscript*. Based on the produced lists of DEGs, the diseases were further analyzed following [IMIDs-analyses.md](./vignettes/IMIDs-analyses.md).

## Environemental set-up to run these scripts

We have used R v4.1.1 and v4.0.4, Python 3.7.6, and MATLAB_R2019b for the development of these codes. 

### Quality assessment and sorting

Package versions:

Seurat 4.0.4

### Clustering and cell type identification by Leiden's algorithm

Package versions:

Scanpy 1.4.5  
DCA 0.2.3  
Tensorflow 1.14.0  
Pandas 0.25.1
Numpy 1.17.2

### Data normalization and Differential expression analysis

Package versions:

scVI 0.7.1  
Pandas 1.2.4  
Numpy 1.21.2  
Scanpy 1.9.1  

### Cell type identification using marker genes

Package versions:

Seurat 4.0.4  
ggplot2 3.3.5  
dplyr 1.0.7  
cowplot 1.1.1  
patchwork 1.1.1  
pheatmap 1.1.12

### Identification of cellular interactions and MCDM/MO-MCDM construction

Package versions:

nichenetr 1.0.0  
dplyr 1.0.7  

### Ranking of URs based on their downstream effect

Package versions:

dplyr 1.0.7  
pheatmap 1.0.12

### Connective Pathway Analysis

Package versions:

pheatmap 1.0.12  
stats 4.0.4

