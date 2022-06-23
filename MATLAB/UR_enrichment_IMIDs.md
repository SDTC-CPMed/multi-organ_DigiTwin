# UR enrichment analysis of IMIDs

```matlab:Code
clear all
close all
clc
warning('off', 'all')
```

# Input requirements
### For this analysis we will need:

   -  The metadata table that summarizes the datasets and where the files for each dataset are stored 
   -  Gene information 
   -  The results from connective pathway analysis 

### We can find the files at following locations:

\hfill \break


```matlab:Code
% DEGs
path_metadata = '../data/Connective_pathway_analysis/AllDatasets_DEGs_martin.mat';
path_gene_info = '../data/Connective_pathway_analysis/gene_info_type_of_gene_2020_03_03.mat';
% connective pathway analysis
path_cpa_program2 = '../data/Connective_pathway_analysis/TreeStructure_nodes2_CLUSTER2_AID_noblood.txt';
path_cpa_program1 = '../data/Connective_pathway_analysis/TreeStructure_nodes2_AID_noblood.txt';
```

### We will produce


A table with combined, FDR corrected pvalues and the count of significant p values across all cell types for each program/subprogram and UR 


### We will save the file here:

\hfill \break


```matlab:Code
path_output = '../data/UR_analysis/UR_predictions_IMIDs_disease_Pvals.xlsx';
```

  
# Prepare the data
  
### Load the metadata table 

\hfill \break


```matlab:Code
%Load the metadata table

load(path_metadata,'AllDatasets')
AllDatasets(strcmp(AllDatasets.active_USE,'unknown'),:)=[]; %WHY DO WE DO THIS//ask Danka
AllDatasets(strcmp(AllDatasets.key,'GSE40568_SS_labial_salivary_gland'),:)=[]; %WHY DO WE DO THIS // ask Danka
head(AllDatasets, 5)
```

| |file_name|file_path|folder_name|folder_path|disease|disease_short|disease_short_broad_categories|tissue|tissue_detail|npatients|ncontrols|gse|active_USE|color_tissue|color_disease_broad|color_disease|deg_key|Disease|disease_key|Datasets|noPatient_Xinxiu|noControl_Xinxiu|Tissues_Xinxiu|Refs|numberOfDEGs|tissues_key|group|key|UR|significantURs|deg_key2|deg_file_name|deg_file_path|DEG|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'UR_GSE16161 skin 9A...|'/Users/danga10/Docu...|'2 Main affected org...|'/Users/danga10/Docu...|'atopic dermatitis'|'AD'|'AD'|'skin'|'skin'|9|9|'GSE16161'|'yes'|'\#eedc1b'|'\#7aff00'|'\#7aff00'|'2 Main affected org...|'Atopic dermatitis'|'AD'|'GSE16161'|9|9|'Skin'|'PMID: 20004782'|27971|'skin'|"2 Main affected org...|'GSE16161_AD_skin'|202x4 table|43x1 cell|'GSE16161_skin_AD'|'GSE16161_Skin AD 9 ...|'/Users/danga10/Docu...|5930x4 table|
|2|'UR_GSE32924_AD_skin...|'/Users/danga10/Docu...|'2 Main affected org...|'/Users/danga10/Docu...|'atopic dermatitis'|'AD'|'AD'|'skin'|'skin'|NaN|8|'GSE32924'|'yes'|'\#eedc1b'|'\#7aff00'|'\#7aff00'|'2 Main affected org...|''|'AD'|'GSE32924'|13|8|'Skin'|'PMID: 21388663'|6939|'skin'|"2 Main affected org...|'GSE32924_AD_skin'|134x4 table|22x1 cell|'GSE32924_skin_AD'|'GSE32924_Skin_AD 13...|'/Users/danga10/Docu...|5997x4 table|
|3|'GSE32924_AD_Inactiv...|'/Users/danga10/Docu...|'4 Non inflamed tiss...|'/Users/danga10/Docu...|'atopic dermatitis'|'AD'|'AD'|'skin_inactive'|'skin_inactive'|NaN|8|'GSE32924'|'no'|'\#eee8a5'|'\#7aff00'|'\#7aff00'|'4 Non inflamed tiss...|'Atopic dermatitis'|'AD'|'GSE32924'|12|8|'Non lesional skin'|'PMID: 21388663'|5715|'skin_inactive'|"4 Non inflamed tiss...|'GSE32924_AD_skin_in...|53x4 table|[{'ESR1',};{'PGR'...|'GSE32924_skin_inact...|'GSE32924_uninflamed...|'/Users/danga10/Docu...|6015x4 table|
|4|'UR_GSE16879_colon_C...|'/Users/danga10/Docu...|'2 Main affected org...|'/Users/danga10/Docu...|'Crohn’s disease'|'CD'|'IBD'|'colon'|'colon'|NaN|6|'GSE16879'|'yes'|'\#c16919'|'\#c16919'|'\#c16919'|'2 Main affected org...|'Crohn’s colitis'|'CD'|'GSE16879'|19|6|'Colonic mucosa'|'PMID: 19956723'|5295|'colon'|"2 Main affected org...|'GSE16879_CD_colon'|225x4 table|96x1 cell|'GSE16879_colon_CD'|'GSE16879_colon_CD.c...|'/Users/danga10/Docu...|5930x4 table|
|5|'UR_GSE179285_ascend...|'/Users/danga10/Docu...|'2 Main affected org...|'/Users/danga10/Docu...|'Crohn’s disease'|'CD'|'IBD'|'colon'|'ascending colon'|35|25|'GSE179285'|'yes'|'\#c16919'|'\#c16919'|'\#c16919'|'2 Main affected org...|''|'CD'|'GSE179285'|14|12|'Ascending/ descendi...|''|7874|'colon'|"2 Main affected org...|'GSE179285_CD_ascend...|192x4 table|92x1 cell|'GSE179285_colon_CD'|'GSE179285_ascending...|'/Users/danga10/Docu...|6012x4 table|

### Load IMID_Ps and IMID_SPs from connective pathway analysis

\hfill \break


```matlab:Code
%% new section
 
load(path_gene_info)
gene_info = unique(gene_info(:,{'GeneID','Symbol','Synonyms','type_of_gene'}));
[Fn, AllDatasets] = preprocess(gene_info, AllDatasets);
clearvars -except Fn savename AllDatasets gene_info path_cpa_program2 path_cpa_program1 path_output
head(Fn, 5)
```

| |key|file_name|group|DEG|broad_label|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'UR_GSE16161 skin 9A...|'UR_GSE16161 skin 9A...|'Inflamed'|16811x1 cell|'AD_active'|
|2|'UR_GSE32924_AD_skin...|'UR_GSE32924_AD_skin...|'Inflamed'|5289x1 cell|'AD_active'|
|3|'GSE32924_AD_Inactiv...|'GSE32924_AD_Inactiv...|'Noninflamed'|4262x1 cell|'AD_inactive'|
|4|'UR_GSE16879_colon_C...|'UR_GSE16879_colon_C...|'Inflamed'|3709x1 cell|'CD_active'|
|5|'UR_GSE179285_ascend...|'UR_GSE179285_ascend...|'Inflamed'|5829x1 cell|'CD_active'|

### Load IMID_SPs and IMID_Ps

\hfill \break


```matlab:Code

%% load SPs %% add subprograms of program 1
SP = load_SPs(path_cpa_program1, path_cpa_program2)
```

| |SP|AllMolecules|
|:--:|:--:|:--:|
|1|'1.1'|1x1108 cell|
|2|'1.10'|1x1116 cell|
|3|'1.2'|1x368  cell|
|4|'1.3'|1x1492 cell|
|5|'1.4'|1x1582 cell|
|6|'1.5'|1x864  cell|
|7|'1.6'|1x2083 cell|
|8|'1.7'|1x1211 cell|
|9|'1.8'|1x1255 cell|
|10|'1.9'|1x305  cell|
|11|'2.1'|1x24   cell|
|12|'2.10'|1x57   cell|
|13|'2.11'|1x79   cell|
|14|'2.12'|1x81   cell|
|15|'2.13'|1x21   cell|
|16|'2.2'|1x1735 cell|
|17|'2.3'|1x229  cell|
|18|'2.4'|1x58   cell|
|19|'2.5'|1x2173 cell|
|20|'2.6'|1x411  cell|
|21|'2.7'|1x298  cell|
|22|'2.8'|1x35   cell|
|23|'2.9'|1x445  cell|
|24|'P1'|1x4693 cell|
|25|'P2'|1x3980 cell|


```matlab:Code
clearvars -except SP Fn savename AllDatasets gene_info path_output

head(SP, 5)
```

| |SP|AllMolecules|
|:--:|:--:|:--:|
|1|'1.1'|1x1108 cell|
|2|'1.10'|1x1116 cell|
|3|'1.2'|1x368  cell|
|4|'1.3'|1x1492 cell|
|5|'1.4'|1x1582 cell|

### Load URs

\hfill \break


```matlab:Code
 
%% load URs: 
UR = load_URs(AllDatasets);
clearvars -except SP Fn savename UR path_output
head(UR, 5)
```

| |UR|key|file_name|group|DS|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'ESR1'|'UR_GSE16161 skin 9A...|'UR_GSE16161 skin 9A...|'Inflamed'|1x735 cell|
|2|'TGFB1'|'UR_GSE16161 skin 9A...|'UR_GSE16161 skin 9A...|'Inflamed'|1x846 cell|
|3|'ESR2'|'UR_GSE16161 skin 9A...|'UR_GSE16161 skin 9A...|'Inflamed'|1x364 cell|
|4|'IFNG'|'UR_GSE16161 skin 9A...|'UR_GSE16161 skin 9A...|'Inflamed'|1x666 cell|
|5|'AGT'|'UR_GSE16161 skin 9A...|'UR_GSE16161 skin 9A...|'Inflamed'|1x415 cell|

  
  
# Main analysis 
### Fisher test enrichment of UR DS in SPs:

\hfill \break


```matlab:Code

FishertestR = Fishers_test(SP,UR,Fn)
```


```text:Output
FishertestR = 
      cellsdes: [25x2 table]
        rowdes: [32x5 table]
    coldes_URs: {506x1 cell}
    moreURinfo: [3129x5 table]
          Pval: {25x1 cell}
            OR: {25x1 cell}

FishertestR = 
      cellsdes: [25x2 table]
        rowdes: [32x5 table]
    coldes_URs: {506x1 cell}
    moreURinfo: [3129x5 table]
          Pval: {25x1 cell}
            OR: {25x1 cell}

```

### combine pvalues over joint and muscle:

\hfill \break


```matlab:Code

CombinedOverInfandNoninf = Fishers_method(FishertestR)
```


```text:Output
CombinedOverInfandNoninf = 
        cols_URs: {506x1 cell}
        rows_SPs: {25x1 cell}
       CombP_Inf: [25x506 double]
    CombP_Noninf: [25x506 double]
       CombP_all: [25x506 double]

```

### FDR correction over Programs

\hfill \break


```matlab:Code

% For each program and UR, do FDR correction 
CombinedOverInfandNoninf = FDR_correction(CombinedOverInfandNoninf)
```


```text:Output
CombinedOverInfandNoninf = 
        cols_URs: {506x1 cell}
        rows_SPs: {25x1 cell}
       CombP_Inf: [25x506 double]
    CombP_Noninf: [25x506 double]
       CombP_all: [25x506 double]
         FDR_Inf: [25x506 double]
      FDR_Noninf: [25x506 double]
         FDR_All: [25x506 double]

```

  
  
### Count significant pvalues over joints and muscles:

\hfill \break


```matlab:Code

% For each program and UR, count how in how many cell types was UR significant 
[count_Inf, count_Noninf] = count_significant(FishertestR);

CombinedOverInfandNoninf.count_Inf = count_Inf;
CombinedOverInfandNoninf.count_Noninf = count_Noninf;

CombinedOverInfandNoninf
```


```text:Output
CombinedOverInfandNoninf = 
        cols_URs: {506x1 cell}
        rows_SPs: {25x1 cell}
       CombP_Inf: [25x506 double]
    CombP_Noninf: [25x506 double]
       CombP_all: [25x506 double]
         FDR_Inf: [25x506 double]
      FDR_Noninf: [25x506 double]
         FDR_All: [25x506 double]
       count_Inf: [25x506 double]
    count_Noninf: [25x506 double]

```

### Combine pvals for each disease

\hfill \break


```matlab:Code
CombinedPval_Disease = combine_Pvals(FishertestR);
```


# Final output

```matlab:Code

CombinedPval  = format_output(CombinedOverInfandNoninf,CombinedPval_Disease)
head(CombinedPval, 5)
```

| |SP|UR|CombinedP_Inflamed|CombinedP_Noninflamed|CombinedP_All|FDR_Inflamed|FDR_Noninflamed|FDR_All|count_Inflamed|count_Noninflamed|AD_active|AD_inactive|CD_active|CD_inactive|JM_active|JM_inactive|PSO_active|PSO_inactive|RA_active|SS_active|SSc_active|UC_active|UC_inactive|at_risk_T1D_inactive|lupus_active|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'1.1'|'ACKR1'|9.9633e-01|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|0|0|1.0000e+00|1|1.8783e-01|1.0000e+00|1|1|3.7841e-01|1|1.0000e+00|1.0000e+00|1.2591e-01|6.7796e-01|1.0000e+00|1|1.0000e+00|
|2|'1.10'|'ACKR1'|3.9770e-03|1.0000e+00|1.5458e-01|1.2736e-02|1.0000e+00|4.5213e-01|7|0|1.0000e+00|1|2.1877e-04|1.0000e+00|1|1|4.5754e-04|1|1.0000e+00|1.0000e+00|6.1443e-04|1.0047e-01|1.0000e+00|1|1.0000e+00|
|3|'1.2'|'ACKR1'|3.5817e-01|1.0000e+00|9.1909e-01|1.0000e+00|1.0000e+00|1.0000e+00|1|0|1.0000e+00|1|3.0478e-02|1.0000e+00|1|1|5.5395e-04|1|1.0000e+00|1.0000e+00|2.6368e-02|4.3447e-01|1.0000e+00|1|1.0000e+00|
|4|'1.3'|'ACKR1'|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|0|0|1.0000e+00|1|1.0000e+00|1.0000e+00|1|1|1.0000e+00|1|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|1|1.0000e+00|
|5|'1.4'|'ACKR1'|9.9971e-01|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|0|0|1.0000e+00|1|3.2724e-01|1.0000e+00|1|1|4.5066e-01|1|1.0000e+00|1.0000e+00|2.0132e-01|7.8095e-01|1.0000e+00|1|1.0000e+00|

  
### Save the data

\hfill \break


```matlab:Code
%writetable(CombinedPval,path_output)

```

