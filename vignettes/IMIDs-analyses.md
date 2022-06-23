# Meta analysis of 11 IMIDs for Connective pathway analysis and UR prioritization

## Process the IPA output (MATLAB code)

**Input requirements**

The enriched pathways produced by IPA including the genes that belong to
these pathways

**We will produce**

The Jaccard Index between the pathways

``` matlab:code
%% Define the main path
InputOutputFiles = '../data/CPA_InputFiles/';

%% Main analyses
% define variables:
savename = 'IMID';
PATHH1 = sprintf('%sPathwayenrichment_results/IMIDs/',InputOutputFiles);
FN = readtable(sprintf('%s/PathFilesDescription_IMIDs.csv',InputOutputFiles));
% read in all pathways, filter significant ones, find all involved genes
% and summarize pathway activation direction:
[AllEnrichedPaths,filesForOverlap] = ReadInIPAPathwayEnrichments(FN,PATHH1);
% calculate jaccard index:
Jaccard = CalculateJaccardIndex(AllEnrichedPaths);

%% Save the outputs:
writetable(Jaccard,sprintf('%sJaccardIndex_%s.csv',InputOutputFiles,savename))
writetable(AllEnrichedPaths,sprintf('%sPathInfo_%s.csv',InputOutputFiles,savename))
writetable(FN,sprintf('%sDatasetInfo_%s.csv',InputOutputFiles,savename),'delimiter','\t')
writetable(filesForOverlap,sprintf('%sDieasesPathways_%s.csv',InputOutputFiles,savename),'delimiter','\t')
```

## Connective pathway analysis (R code)

Similarly as described for the [CIA analysis - Connective Pathway
Analysis](./CIA-analyses.md), we have performed connective pathway
analysis for the IMID datasets.

``` r
source("../R/CPA_UR_rankings_functions.R")

MainPath = "data/CPA_InputFiles/"

# load pre-computed Jaccard index matrix and pathway information plus define name which is to be used to save results:

xx = as.matrix(read.csv(paste(MainPath, 'JaccardIndex_IMID.csv',sep='')))
pathinfo = read.table(paste(MainPath,'PathInfo_IMID.csv',sep=''),sep=',',header=1) 
SAVENAME = 'IMID'

pathinfo = runCPA(xx,pathinfo,SAVENAME,1.78)
# save results to file:
write.table(file=paste(MainPath,'CPA_',SAVENAME,'_all_pathways.txt',sep=''),pathinfo,sep="\t", col.names=T, row.names=F, quote=F)

# count the ratios of pathways showing same and opposing activation pattern in inflamed vs non-inflamed organ:
ratios = count.same.and.opposing.activations(pathinfo)
```

![Fig_tree_pathways_AID](IMIDs-analyses_files/figure-markdown_github/Fig_tree_pathways_AID.png)

## Overlap of the CPA of IMIDs vs individual IMIDs (R code)

In order to get a better overview on which programs and sub-programs of
CPA are enriched in individual diseases, we calulated a Fisher Exact
Test. To test if the programs (IMID_P1 and IMID_P2) derived from all
analyzed IMIDs overlapped with programs from individual IMIDs in
inflamed and non-inflamed organ sites, separately, we performed Fisher’s
exact tests (right tailed), using all pathways in connective pathway
analysis as a background, followed by correction for multiple testing
using the Benjamini-Hochberg procedure. These analyses were repeated for
IMID_subprograms (IMID_SPs) and subprograms from individual IMIDs.
Disease pathways were defined as all pathways significantly enriched in
a particular disease and inflammation state (IPA, p \< 0.05). In case
two or more datasets were representative of the same disease and
condition (for example UC inflamed organs), pathway enrichment p values
were combined with Fisher’s method. Enriched pathways were considered
those whose combined p \< 0.05

``` r
# prepare files to check overlap between programs and enriched pathways per IMID:
Dis.Pval = read.table(paste(MainPath,'DieasesPathways_IMID.csv',sep=''),sep='\t',quote = "",header=T)
OverlapingPaths = FindOverlapingPaths(Dis.Pval,pathinfo)

source("../R/plot_overlap.R")
overlap = read.csv('../data/IMIDs_pathway_overlap_with_SPs_reshaped_for_dot_plotALL.txt')

temp_plot <- plot_overlap(overlap)
temp_plot
```

![](IMIDs-analyses_files/figure-markdown_github/unnamed-chunk-3-1.png)

## UR enrichment analysis of IMIDs (MATLAB code)

**Input requirements**

-   The metadata table that summarizes the datasets and where the files
    for each dataset are stored
-   Gene information
-   The results from connective pathway analysis

**We will produce**

A table with combined, FDR corrected pvalues and the count of
significant p values across all cell types for each program/subprogram
and UR

``` matlab:code
%% INPUT
% DEGs
path_metadata = '../data/Connective_pathway_analysis/AllDatasets_DEGs_martin.mat';
path_gene_info = '../data/Connective_pathway_analysis/gene_info_type_of_gene_2020_03_03.mat';
% connective pathway analysis
path_cpa_program2 = '../data/Connective_pathway_analysis/TreeStructure_nodes2_CLUSTER2_AID_noblood.txt';
path_cpa_program1 = '../data/Connective_pathway_analysis/TreeStructure_nodes2_AID_noblood.txt';

%% OUTPUT_path
path_output = '../data/UR_analysis/UR_predictions_IMIDs_disease_Pvals.xlsx';

%%Load the metadata table
load(path_metadata,'AllDatasets')

%% Load IMID_Ps and IMID_SPs from connective pathway analysis
load(path_gene_info)
gene_info = unique(gene_info(:,{'GeneID','Symbol','Synonyms','type_of_gene'}));
[Fn, AllDatasets] = preprocess(gene_info, AllDatasets);

%% load SPs %% add subprograms of program 1
SP = load_SPs(path_cpa_program1, path_cpa_program2)

%% Load URs: 
UR = load_URs(AllDatasets);
clearvars -except SP Fn savename UR path_output

%% Fisher test enrichment of UR DS in SPs:
FishertestR = Fishers_test(SP,UR,Fn)

%% Combine pvalues over joint and muscle:
CombinedOverInfandNoninf = Fishers_method(FishertestR)

%% For each program and UR, do FDR correction 
CombinedOverInfandNoninf = FDR_correction(CombinedOverInfandNoninf)

%% For each program and UR, count how in how many cell types was UR significant 
[count_Inf, count_Noninf] = count_significant(FishertestR);
CombinedOverInfandNoninf.count_Inf = count_Inf;
CombinedOverInfandNoninf.count_Noninf = count_Noninf;

%% Combine pvals for each disease
CombinedPval_Disease = combine_Pvals(FishertestR);

%% Format output
CombinedPval  = format_output(CombinedOverInfandNoninf,CombinedPval_Disease)

%% Save the data
writetable(CombinedPval,path_output)
head(CombinedPval, 5)
```

|     |   SP   |   UR    | CombinedP_Inflamed | CombinedP_Noninflamed | CombinedP_All | FDR_Inflamed | FDR_Noninflamed |  FDR_All   | count_Inflamed | count_Noninflamed | AD_active  | AD_inactive | CD_active  | CD_inactive | JM_active | JM_inactive | PSO_active | PSO_inactive | RA_active  | SS_active  | SSc_active | UC_active  | UC_inactive | at_risk_T1D_inactive | lupus_active |
|:---:|:------:|:-------:|:------------------:|:---------------------:|:-------------:|:------------:|:---------------:|:----------:|:--------------:|:-----------------:|:----------:|:-----------:|:----------:|:-----------:|:---------:|:-----------:|:----------:|:------------:|:----------:|:----------:|:----------:|:----------:|:-----------:|:--------------------:|:------------:|
|  1  | ‘1.1’  | ‘ACKR1’ |     9.9633e-01     |      1.0000e+00       |  1.0000e+00   |  1.0000e+00  |   1.0000e+00    | 1.0000e+00 |       0        |         0         | 1.0000e+00 |      1      | 1.8783e-01 | 1.0000e+00  |     1     |      1      | 3.7841e-01 |      1       | 1.0000e+00 | 1.0000e+00 | 1.2591e-01 | 6.7796e-01 | 1.0000e+00  |          1           |  1.0000e+00  |
|  2  | ‘1.10’ | ‘ACKR1’ |     3.9770e-03     |      1.0000e+00       |  1.5458e-01   |  1.2736e-02  |   1.0000e+00    | 4.5213e-01 |       7        |         0         | 1.0000e+00 |      1      | 2.1877e-04 | 1.0000e+00  |     1     |      1      | 4.5754e-04 |      1       | 1.0000e+00 | 1.0000e+00 | 6.1443e-04 | 1.0047e-01 | 1.0000e+00  |          1           |  1.0000e+00  |
|  3  | ‘1.2’  | ‘ACKR1’ |     3.5817e-01     |      1.0000e+00       |  9.1909e-01   |  1.0000e+00  |   1.0000e+00    | 1.0000e+00 |       1        |         0         | 1.0000e+00 |      1      | 3.0478e-02 | 1.0000e+00  |     1     |      1      | 5.5395e-04 |      1       | 1.0000e+00 | 1.0000e+00 | 2.6368e-02 | 4.3447e-01 | 1.0000e+00  |          1           |  1.0000e+00  |
|  4  | ‘1.3’  | ‘ACKR1’ |     1.0000e+00     |      1.0000e+00       |  1.0000e+00   |  1.0000e+00  |   1.0000e+00    | 1.0000e+00 |       0        |         0         | 1.0000e+00 |      1      | 1.0000e+00 | 1.0000e+00  |     1     |      1      | 1.0000e+00 |      1       | 1.0000e+00 | 1.0000e+00 | 1.0000e+00 | 1.0000e+00 | 1.0000e+00  |          1           |  1.0000e+00  |
|  5  | ‘1.4’  | ‘ACKR1’ |     9.9971e-01     |      1.0000e+00       |  1.0000e+00   |  1.0000e+00  |   1.0000e+00    | 1.0000e+00 |       0        |         0         | 1.0000e+00 |      1      | 3.2724e-01 | 1.0000e+00  |     1     |      1      | 4.5066e-01 |      1       | 1.0000e+00 | 1.0000e+00 | 2.0132e-01 | 7.8095e-01 | 1.0000e+00  |          1           |  1.0000e+00  |

## logFC and z-score analyses (Python code)

**Input requirements**

-   The UR enrichment scores computed above
-   Folder including differentially expressed genes for each dataset
-   Folder including z scores for each IPA result

**We will produce**

Dot plots showing the A) predicted activities and B) fold changes, of
the shared URs of IMID_SP1.6. (Figure S8)

``` python
## Input
URs_all_diseases = pd.read_table('../data/UR_analysis/UR_predictions_IMIDs_disease_Pvals.txt', sep = ',')
path_DEGs = '../data/AllDEGfilesMovedToOneFolder/'
path_z_scores = '../data/UR_analysis/z_scores/ '

## Output paths
path_Data_S15 = '../data/UR_analysis/Data S15.xlsx'
path_URs_logFC = '../data/UR_analysis/UR_IMID_summary_logFC.csv'
path_URs_zScore = '../data/UR_analysis/UR_IMID_summary_z.csv'

## Preprocess the data
URs_all_diseases = preprocess_the_input(URs_all_diseases)
Datasets_Inf, Datasets_Noninf = set_up_structure_for_output_files()

## logFC_analysis
summary_logFC = logFC_analysis(Datasets_Inf, Datasets_Noninf, path_DEGs)
summary_logFC.to_csv(path_URs_logFC)

## z score analysis
summary_zScores = zScore_analysis(Datasets_Inf, Datasets_Noninf, path_z_scores)
summary_zScores.to_csv(path_URs_zScores)
```

## Plot the logFC and z_scores of URs (R code)

**URs z score**

``` r
source('../R/plot_zScore.R')
URs_zScore = read.csv('../data/UR_analysis/UR_IMID_summary_z.csv')
temp_plot <- plot_zScore(URs_zScore)
temp_plot
```

![](IMIDs-analyses_files/figure-markdown_github/unnamed-chunk-6-1.png)

**URs logFC**

``` r
source("../R/plot_logFC.R")
URs_logFC = read.csv('../data/UR_analysis/UR_IMID_summary_logFC.csv')
temp_plot <- plot_logFC(URs_logFC)
temp_plot
```

![](IMIDs-analyses_files/figure-markdown_github/unnamed-chunk-7-1.png)
